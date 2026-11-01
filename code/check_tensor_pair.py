"""Analytic polynomial controls for non-vacuum pair cubature and bounds."""
from pathlib import Path
import sys,json,hashlib
root=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(root/'tmp/arb-runtime')]
import numpy as np
from flint import arb,acb,ctx
from tensor_pair_overlap import PairOverlap,fourier
from chirp_residual_certificate import complex_exact,interval,norm_squared
from residual_overlap_rectangle import phase_product

def moments(a,b,count):
    a,b=arb(a),arb(b);ea=(-a*a).exp();eb=(-b*b).exp();p=arb.pi().sqrt()
    m=[(b.erf()-a.erf())/2,(ea-eb)/(2*p)]
    for k in range(2,count):m.append((a**(k-1)*ea-b**(k-1)*eb)/(2*p)+(arb(k-1)/2)*m[k-2])
    return m

def main():
    ctx.prec=320;checks=[]
    def record(name,ok,**details):checks.append(dict(name=name,passed=bool(ok),**details))
    v=np.zeros((4,4),complex);w=v.copy()
    v[1,2]=.5+.25j;v[3,0]=-.75j;v[0,1]=.25
    w[2,3]=.125-.5j;w[1,0]=.5;w[0,1]=1+.5j
    ev=PairOverlap(v,w,terms=3,order=24,unit=True)
    ps=[[arb(1)],[arb(0),arb(2).sqrt()],[-1/arb(2).sqrt(),arb(0),arb(2).sqrt()],
        [arb(0),-12/arb(48).sqrt(),arb(0),8/arb(48).sqrt()]]
    def gram(a,b):
        m=moments(a,b,7)
        return [[sum((ps[i][k]*ps[j][l]*m[k+l] for k in range(len(ps[i])) for l in range(len(ps[j]))),arb(0)) for j in range(4)] for i in range(4)]
    for rect in [(-1,1,-1,1),(-.5,1.25,.25,1.5),(-2,-.25,-1,.5)]:
        gx,gy=gram(*rect[:2]),gram(*rect[2:]);exact=acb(0)
        for i,j in np.argwhere(w!=0):
            for k,l in np.argwhere(v!=0):exact+=complex_exact(w[i,j]).conjugate()*complex_exact(v[k,l])*gx[i][k]*gy[j][l]
        c=ev.bound(rect);q=ev.quadrature(rect)
        record('nonvacuum antiderivative',c['kind']=='gauss' and abs(q-exact)<c['error'],rect=rect,error=interval(c['error']),difference=interval(abs(q-exact)))
    for c in [v,w]:record('Fourier preserves exact dyadic norm',(norm_squared(c)-norm_squared(fourier(c))).is_zero())
    ev=PairOverlap(v,w,terms=8,order=24)
    for x,y,b in [(.25,.5,.125),(3,2,-.125),(-3,2,.25),(10,-4,.125)]:
        z=acb(arb(x),arb(b));t=acb(arb(y))
        exact=abs(phase_product(arb.pi().sqrt()*z,arb.pi().sqrt()*t,8))
        bound=ev.phase_bound(z,t)
        record('phase inequality at complex point',exact<bound,point=[x,y,b])
    pole=arb.pi().sqrt()
    record('pole rejection',not ev.phase_bound(acb(pole),acb(pole)).is_finite())
    pure=np.zeros((4,4),complex);pure[3,2]=1
    pv=PairOverlap(pure,pure,unit=True);L=arb(4)
    m=[L.erfc()]
    for k in range(1,4):m.append(L**(2*k-1)*(-L*L).exp()/arb.pi().sqrt()+arb(2*k-1)/2*m[-1])
    tx=arb(4)/3*m[3]-4*m[2]+3*m[1];ty=2*m[2]-2*m[1]+m[0]/2
    exact_tail=tx+ty-tx*ty
    record('nonvacuum Hermite tail via exact recurrence',exact_tail>0 and exact_tail<pv.tails(4)[0],tail=interval(exact_tail),bound=interval(pv.tails(4)[0]))
    result=dict(passed=all(c['passed'] for c in checks),checks=checks,
        scope='Controls for explicit general pair formulas; no complete F chain certified here.',
        hashes={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in ['code/check_tensor_pair.py','code/tensor_pair_overlap.py']})
    with (root/'results/stage4/residual-certificate/tensor-pair-controls.json').open('x') as f:json.dump(result,f,indent=2)
    print(json.dumps(dict(passed=result['passed'],count=len(checks),failed=[c for c in checks if not c['passed']])))
    assert result['passed']

if __name__=='__main__':main()
