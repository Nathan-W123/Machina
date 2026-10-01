"""Extract vector polylines (with stroke colour and CTM) and text positions from a PDF page.
Usage: python3 extract_paths.py pdf page_index(0-based) out.json"""
import sys, json
sys.modules['cryptography'] = None
import pypdf
from pypdf.generic import ContentStream

def mul(a, b):  # 2D affine [a b c d e f]
    return [a[0]*b[0]+a[1]*b[2], a[0]*b[1]+a[1]*b[3], a[2]*b[0]+a[3]*b[2], a[2]*b[1]+a[3]*b[3],
            a[4]*b[0]+a[5]*b[2]+b[4], a[4]*b[1]+a[5]*b[3]+b[5]]
def ap(m, x, y): return (m[0]*x+m[2]*y+m[4], m[1]*x+m[3]*y+m[5])

pdf, pg, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
r = pypdf.PdfReader(pdf); page = r.pages[pg]
cs = ContentStream(page.get_contents(), r)
ctm = [1,0,0,1,0,0]; stack = []; stroke = None; fill = None; lw = 1
paths = []; cur = []; sub = []
texts = []; tm = [1,0,0,1,0,0]; tlm = [1,0,0,1,0,0]; lead = 0
for ops, op in cs.operations:
    op = op.decode() if isinstance(op, bytes) else op
    f = [float(x) for x in ops] if op in ('cm','m','l','c','v','y','re','RG','rg','K','k','G','g','w','Tm','Td','TD','TL','SC','SCN','sc','scn') and all(hasattr(x,'__float__') for x in ops) else None
    if op == 'q': stack.append((ctm[:], stroke, fill, lw))
    elif op == 'Q' and stack: ctm, stroke, fill, lw = stack.pop()
    elif op == 'cm' and f: ctm = mul(f, ctm)
    elif op in ('RG','K','G','SC','SCN') and f is not None: stroke = tuple(round(v,3) for v in f)
    elif op in ('rg','k','g','sc','scn') and f is not None: fill = tuple(round(v,3) for v in f)
    elif op == 'w' and f: lw = f[0]
    elif op == 'm' and f:
        if len(sub) > 1: cur.append(sub)
        sub = [ap(ctm, f[0], f[1])]
    elif op == 'l' and f: sub.append(ap(ctm, f[0], f[1]))
    elif op in ('c',) and f: sub.append(ap(ctm, f[4], f[5]))
    elif op in ('v','y') and f: sub.append(ap(ctm, f[2], f[3]))
    elif op == 're' and f:
        x,y,w,h = f; 
        if len(sub) > 1: cur.append(sub)
        cur.append([ap(ctm,x,y),ap(ctm,x+w,y),ap(ctm,x+w,y+h),ap(ctm,x,y+h),ap(ctm,x,y)]); sub=[]
    elif op in ('h',): pass
    elif op in ('S','s','f','F','f*','B','B*','b','b*','n'):
        if len(sub) > 1: cur.append(sub)
        if op != 'n' and cur:
            paths.append({'op': op, 'stroke': stroke, 'fill': fill, 'lw': lw*abs(ctm[0]), 'subs': cur})
        cur = []; sub = []
    elif op == 'BT': tm = [1,0,0,1,0,0]; tlm = tm[:]
    elif op == 'Tm' and f: tm = f[:]; tlm = f[:]
    elif op in ('Td','TD') and f:
        tlm = mul([1,0,0,1,f[0],f[1]], tlm); tm = tlm[:]
        if op == 'TD': lead = -f[1]
    elif op == 'T*': tlm = mul([1,0,0,1,0,-lead], tlm); tm = tlm[:]
    elif op in ('Tj','TJ',"'"):
        a = ops[0]
        s = ''.join(x for x in a if isinstance(x, str)) if op == 'TJ' else str(a)
        x, y = ap(mul(tm, ctm), 0, 0)
        texts.append({'x': x, 'y': y, 's': s})
json.dump({'paths': paths, 'texts': texts, 'mediabox': [float(v) for v in page.mediabox]}, open(out, 'w'))
print(len(paths), 'paths', len(texts), 'texts')
