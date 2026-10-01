"""Build two challenge decks on the convergence lens' 0.833 mm IM5 TC deck
(cases/h0.833_L1_im_tp5, which is the base for comparison):
  kin_h0833 : + material lens' kin_shutov material (2 AF backstresses)
  dsif_h0833: + process lens' DSIF support tool (dsif_ng support.csv)"""
import json, shutil
from pathlib import Path
A = Path(__file__).resolve().parent.parent
base = A / "convergence/cases/h0.833_L1_im_tp5/deck"
kin_deck = json.loads((A / "material/runs/truncated_cone-s2026-0000/kin_shutov/deck.json").read_text())
ds_dir = A / "process/runs/truncated_cone-s2026-0000/dsif_ng"
ds_deck = json.loads((ds_dir / "deck.json").read_text())
for name in ["kin_h0833", "dsif_h0833"]:
    d = A / "challenge/cases" / name / "deck"
    if d.exists():
        continue
    shutil.copytree(base, d)
    deck = json.loads((d / "deck.json").read_text())
    if name == "kin_h0833":
        deck["material"] = kin_deck["material"]
    else:
        deck["forming"]["tools"] = ds_deck["forming"]["tools"]
        deck["forming"]["steps"][0]["tools"] = ds_deck["forming"]["steps"][0]["tools"]
        shutil.copy(ds_dir / "support.csv", d / "support.csv")
    (d / "deck.json").write_text(json.dumps(deck, indent=1))
    print(name, "written")
