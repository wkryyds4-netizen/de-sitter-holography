"""Pilot: a finite-product overlap on one pole-free rectangle, using Arb.

This is NOT the full-plane overlap. The outer analyticity guard is essential.
No absolute value or complex conjugation of variables is used in callbacks.
"""
from pathlib import Path
import sys,json,time,hashlib,platform
root=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(root/'tmp/arb-runtime'))
from flint import arb,acb,acb_poly,ctx
import numpy as np
import flint
from chirp_residual_certificate import polynomial_matrix,interval

def phase_product(x,y,terms):
    val=acb(1);ii=acb(0,1)
    for j in range(terms):
        a=(2*j+1)*arb.pi()
        val*=(1+(x-a-ii*y).exp())/(1+(x-a+ii*y).exp())
        val*=(1+(y-a-ii*x).exp())/(1+(y-a+ii*x).exp())
    return val

def outer_guard(y):
    # For real inner coordinate x, zeros of either product denominator
    # require Re(y) to be an odd multiple of sqrt(pi). Holomorphy follows
    # uniformly on the compact x interval when all these rows are avoided.
    t=y.real/arb.pi().sqrt()
    if not t.is_finite():return False
    lo,hi=t.lower().floor(),t.upper().ceil()
    if hi-lo>1000:return False
    if not lo.is_exact() or not hi.is_exact():return False
    return not any(t.contains(k) for k in range(int(lo.fmpq().p),int(hi.fmpq().p)+1) if k%2)

def integrate_rectangle(v,w,xa,xb,ya,yb,terms=8):
    a,b=polynomial_matrix(v),polynomial_matrix(w)
    n=v.shape[0]
    pa=[acb_poly([a[i,j] for j in range(n)]) for i in range(n)]
    pb=[acb_poly([b[i,j].conjugate() for j in range(n)]) for i in range(n)]
    pi=arb.pi();scale=pi.sqrt();calls=dict(outer=0,inner=0,rejected_outer=0)
    def outer(y,analytic):
        calls['outer']+=1
        if not outer_guard(y):
            calls['rejected_outer']+=1
            return acb('nan')
        pv=acb_poly([p(y) for p in pa]);pw=acb_poly([p(y) for p in pb])
        def inner(x,analytic_x):
            calls['inner']+=1
            return pw(x)*pv(x)*(-x*x-y*y).exp()*phase_product(scale*x,scale*y,terms)/pi
        return acb.integral(inner,xa,xb,abs_tol=arb('1e-13'),rel_tol=arb('1e-13'),
                            eval_limit=10000,depth_limit=60)
    value=acb.integral(outer,ya,yb,abs_tol=arb('1e-9'),rel_tol=arb('1e-9'),
                       eval_limit=10000,depth_limit=60)
    return value,calls

def main():
    out=root/'results/stage4/residual-certificate/rectangle-pilot.json'
    if out.exists():raise RuntimeError('Do not overwrite')
    ctx.prec=512
    source=root/'results/stage4/flip-domain/coefficients-03.npz'
    with np.load(source,allow_pickle=False) as data:
        v,w=data['right'][0:2]
    start=time.perf_counter()
    value,calls=integrate_rectangle(v,w,-1,1,-1,1)
    assert value.is_finite()
    # Only the compact square [-1,1]^2 in z=x/sqrt(pi), t=y/sqrt(pi).
    # The infinite-product error on this square is bounded uniformly by tail.
    r=(arb.pi().sqrt()-17*arb.pi()).exp()
    tail=4*r/((1-(-2*arb.pi()).exp())*(1-r))
    from chirp_residual_certificate import norm_squared
    overlap_tail=tail*(norm_squared(v)*norm_squared(w)).sqrt()
    data=dict(status='one compact rectangle enclosed; full-plane overlap uncomputed',
        source_file=source.relative_to(root).as_posix(),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        branch='right',step=1,n=64,precision_bits=512,dimensionless_rectangle=[-1,1,-1,1],
        finite_product_terms=8,finite_product_overlap_real=interval(value.real),
        finite_product_overlap_imag=interval(value.imag),infinite_product_error_bound=interval(overlap_tail),
        calls=calls,seconds=time.perf_counter()-start,
        versions=dict(python=platform.python_version(),python_flint=flint.__version__),
        hashes={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in
                ['code/residual_overlap_rectangle.py','code/chirp_residual_certificate.py']},
        missing='Complement of the square, singularity strips, tails, and the remaining F-step overlaps. This pilot is not a continuum pentagon certificate.')
    out.write_text(json.dumps(data,indent=2))
    print(json.dumps(data),flush=True)

if __name__=='__main__':main()
