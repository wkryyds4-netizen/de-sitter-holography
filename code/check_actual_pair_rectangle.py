"""Actual n128 pair: polynomial/nested-Arb versus Hermite/tensor GL."""
from pathlib import Path
import sys,json,hashlib,time
root=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(root/'tmp/arb-runtime'))
import numpy as np
from flint import arb,acb,ctx
from tensor_pair_overlap import PairOverlap,fourier
from residual_overlap_rectangle import integrate_rectangle
from chirp_residual_certificate import interval

def main():
    out=root/'results/stage4/residual-certificate/actual-pair-rectangle.json'
    if out.exists():raise RuntimeError('Preserve old output')
    source=root/'results/stage4/flip-domain/coefficients-06.npz'
    with np.load(source,allow_pickle=False) as data:v,w=data['right'][4:6]
    v,w=fourier(v),fourier(w);ctx.prec=512;start=time.perf_counter()
    # This square avoids every odd-sqrt(pi) pole row required by the nested guard.
    rect=(-.5,.5,-.5,.5)
    nested,calls=integrate_rectangle(v,w,*rect,terms=8)
    ev=PairOverlap(v,w,terms=8,order=24);cells,error=ev.partition(.5,arb('1e-10'))
    tensor=sum((ev.quadrature(c['rect']) if c['kind']=='gauss' else acb(0) for c in cells),acb(0))
    checks=[]
    for part in ['real','imag']:
        a,b=getattr(tensor,part),getattr(nested,part)
        checks.append(dict(part=part,passed=bool(a.lower()-error<b.upper() and a.upper()+error>b.lower())))
    payload=dict(passed=all(x['passed'] for x in checks),branch='right',step=5,n=128,rectangle=rect,terms=8,precision_bits=ctx.prec,
        nested_real=interval(nested.real),nested_imag=interval(nested.imag),tensor_real=interval(tensor.real),tensor_imag=interval(tensor.imag),tensor_error=interval(error),checks=checks,calls=calls,seconds=time.perf_counter()-start,
        scope='Local real-data cross-check only; full-plane proof uses its separate tail and cell certificate.',
        hashes={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in ['code/check_actual_pair_rectangle.py','code/tensor_pair_overlap.py','code/residual_overlap_rectangle.py',source.relative_to(root).as_posix()]})
    with out.open('x') as f:json.dump(payload,f,indent=2)
    print(json.dumps(dict(passed=payload['passed'],seconds=payload['seconds'],calls=calls,nested=nested.str(12),tensor=tensor.str(12))),flush=True)

if __name__=='__main__':main()
