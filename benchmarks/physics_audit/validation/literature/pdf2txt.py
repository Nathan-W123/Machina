"""Extract text from PDFs with pypdf (the system 'cryptography' package panics on import; block it)."""
import sys
sys.modules['cryptography'] = None
import pypdf

for p in sys.argv[1:]:
    r = pypdf.PdfReader(p)
    out = p.rsplit('.', 1)[0] + '.txt'
    with open(out, 'w') as f:
        for i, page in enumerate(r.pages):
            f.write(f"\n\n===== page {i+1} =====\n")
            f.write(page.extract_text() or '')
    print(out, len(r.pages))
