#!/usr/bin/env bash
# Run the stage n26 ML benchmark on a local machine (Ubuntu, or WSL2 Ubuntu on
# Windows). One command, safe to re-run: every step resumes where it stopped.
#
#   bash benchmarks/springback_fine_ml/laptop.sh            # all cores sensible default
#   WORKERS=8 bash benchmarks/springback_fine_ml/laptop.sh  # simulations at a time
#
# Steps: system packages -> Python venv -> build sparlab_form -> unpack the
# stage n18 data set (laptop_seed.tgz: 128 SparLab samples + the 16 test-part
# baseline runs, so nothing already simulated is simulated again) -> run
# stage n26 -> commit the results (you push).
set -euo pipefail
cd "$(dirname "$0")/../.."
REPO=$PWD
B=benchmarks/springback_fine_ml
W=$B/work

cores=$(nproc)
WORKERS=${WORKERS:-$(( cores / 3 > 1 ? cores / 3 : 1 ))}
THREADS=${THREADS:-4}
echo "== $cores logical cores, $WORKERS simulations at a time"

if ! dpkg -s libsuitesparse-dev >/dev/null 2>&1; then
  echo "== installing system packages (asks for your password once)"
  sudo apt-get update -qq
  sudo apt-get install -y build-essential cmake ninja-build libeigen3-dev \
       libsuitesparse-dev libopenblas-dev python3-venv python3-dev git
fi

if [ ! -x .venv/bin/python ]; then
  echo "== Python environment"
  python3 -m venv .venv
  .venv/bin/pip install -q --upgrade pip
  .venv/bin/pip install -q torch --index-url https://download.pytorch.org/whl/cpu
  .venv/bin/pip install -q -e .
fi

if [ ! -x build/bin/sparlab_form ]; then
  echo "== building sparlab_form"
  cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DSPARLAB_BUILD_TESTS=OFF
  cmake --build build --target sparlab_form -j "$cores"
fi
build/bin/sparlab_form --version

if [ ! -f "$W/data/index.jsonl" ]; then
  echo "== unpacking the stage n18 data set"
  mkdir -p "$W"
  tar -xzf "$B/laptop_seed.tgz" -C "$W"
fi

echo "== running stage n26 (log: $W/ml.log); takes hours - keep the laptop plugged in"
echo "   and sleep off. Interrupted? Re-run this script; it resumes."
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH="$REPO/python" \
  .venv/bin/python python/scripts/springback_fine_ml.py --out "$B" --work "$W" \
  --stages 26 --workers "$WORKERS" --threads "$THREADS" 2>&1 | tee -a "$W/ml.log"

if [ -d "$B/stage_n26" ]; then
  git add "$B/stage_n26" "$B/data_status.csv" "$B/config.json" 2>/dev/null || true
  git -c user.name="${GIT_AUTHOR_NAME:-$(git config user.name || echo laptop)}" \
      commit -q -m "Add stage n26 of the fine-mesh ML benchmark (run locally)" || true
  echo "== done. Results: $B/stage_n26/  ->  now run: git push"
fi
