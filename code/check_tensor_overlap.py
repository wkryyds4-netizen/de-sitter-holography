"""Independent low-Hermite antiderivatives and earlier nested-Arb control."""
from pathlib import Path
import sys,json,hashlib
root=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(root/'tmp/arb-runtime'))
import numpy as np
from flint import arb,acb,ctx,fmpq
from tensor_overlap import GaussianOverlap
from chirp_residual_certificate import complex_exact,interval

def moments(a,b):
    a,b=arb(a),arb(b);ea=(-a*a).exp();eb=(-b*b).exp();pi=arb.pi()
    return [(b.erf()-a.erf())/2,(ea-eb)/(2*pi).sqrt(),
            (a*ea-b*eb)/(2*pi).sqrt(),
            ((4*a*a-2)*ea-(4*b*b-2)*eb)/(48*pi).sqrt()]

def main():
    output=root/'results/stage4/residual-certificate/tensor-controls.json'
    if output.exists():raise RuntimeError('Preserve existing controls')
    ctx.prec=256;rows=[]
    w=np.zeros((4,4),complex)
    w[0,0]=1+.5j;w[1,2]=-.25+.75j;w[3,1]=.125-.5j;w[2,3]=.5+.25j
    ev=GaussianOverlap(w,unit=True)
    for rect in [(-1,1,-1,1),(-.5,1.25,.25,1.5),(-2,-.25,-1,.5)]:
        x,y=moments(*rect[:2]),moments(*rect[2:])
        exact=sum((complex_exact(w[i,j]).conjugate()*x[i]*y[j] for i in range(4) for j in range(4)),acb(0))
        cell=ev.bound(rect);q=ev.quadrature(rect)
        ok=cell['kind']=='gauss' and bool(abs(q-exact)<cell['error'])
        rows.append(dict(control='complex nonsymmetric low Hermite',rect=rect,passed=ok,
                         discrepancy=interval(abs(q-exact)),error=interval(cell['error'])))
    src=root/'results/stage4/flip-domain/coefficients-03.npz'
    with np.load(src,allow_pickle=False) as data:w=data['right'][1]
    ev=GaussianOverlap(w,terms=8,order=24)
    # Separate cells give useful analytic bounds at n=64.
    cells,error=ev.partition(1,arb('1e-10'))
    q=sum((ev.quadrature(c['rect']) if c['kind']=='gauss' else acb(0) for c in cells),acb(0))
    pilot_path=root/'results/stage4/residual-certificate/rectangle-pilot.json'
    pilot=json.loads(pilot_path.read_text())
    for part in ['real','imag']:
        item=pilot['finite_product_overlap_'+part]
        lo,hi=[arb(fmpq(int(a),int(b))) for a,b in [item['lower'],item['upper']]]
        result=getattr(q,part)
        # Stronger than interval overlap: new enclosure contains the old one.
        ok=bool(result.upper()-error<lo and result.lower()+error>hi)
        rows.append(dict(control='contains nested Arb pilot '+part,passed=ok,
                         value=interval(result),analytic_error=interval(error),leaves=len(cells)))
    payload=dict(passed=all(r['passed'] for r in rows),checks=rows,precision_bits=ctx.prec,
                 scope='Implementation controls; not proof of all five full-plane F overlaps',
                 hashes={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in
                         ['code/check_tensor_overlap.py','code/tensor_overlap.py',
                          'results/stage4/residual-certificate/rectangle-pilot.json',src.relative_to(root).as_posix()]})
    output.write_text(json.dumps(payload,indent=2))
    print(json.dumps(payload),flush=True)

if __name__=='__main__':main()
