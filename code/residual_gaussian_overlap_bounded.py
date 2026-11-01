"""One full F overlap, exploiting an exactly Gaussian input.

Uses existing Arb integration and an explicit outer analyticity guard.
The other four F overlaps are not computed here.
"""
from pathlib import Path
import sys,json,hashlib,time,platform
root=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(root/'tmp/arb-runtime'))
import numpy as np
import flint
from flint import arb,acb,acb_poly,ctx
from chirp_residual_certificate import polynomial_matrix,norm_squared,interval
from residual_overlap_rectangle import phase_product,outer_guard
from gaussian_outer_bound import gaussian_outer_bound

def main():
    base=root/'results/stage4/residual-certificate'
    out=base/'gaussian-bounded.json';log=base/'gaussian-bounded-segments.jsonl'
    if out.exists() or log.exists():raise RuntimeError('Do not overwrite or restart')
    ctx.prec=768;L=6;rho=arb('1e-8');terms=8
    source=root/'results/stage4/flip-domain/coefficients-06.npz'
    with np.load(source,allow_pickle=False) as data:v,w=data['right'][0:2]
    vacuum=np.zeros_like(v);vacuum[0,0]=1
    assert np.array_equal(v,vacuum)
    nw=norm_squared(w);normw=nw.sqrt();pi=arb.pi();scale=pi.sqrt()
    # The union of two Gaussian coordinate tails has mass <= 2 erfc(L).
    tail=(2*arb(L).erfc()*nw).sqrt()
    r=arb(2*w.shape[0]-1)/(2*w.shape[0])
    kernel=r**(-(w.shape[0]-1))/(pi*(1-r*r)).sqrt()
    poles=[arb(k)*scale for k in [-3,-1,1,3]]
    strips=2*rho*len(poles)*kernel*normw
    first=(L*scale-(2*terms+1)*pi).exp()
    product=4*first/((1-(-2*pi).exp())*(1-first))*normw
    error=tail+strips+product
    b=polynomial_matrix(w);n=w.shape[0]
    pb=[acb_poly([b[i,j].conjugate() for j in range(n)]) for i in range(n)]
    intervals=[];a=arb(-L)
    for p in poles:
        intervals.append((a,p-rho));a=p+rho
    intervals.append((a,arb(L)))
    metadata=dict(status='running',index=6,branch='right',step=1,n=128,precision_bits=ctx.prec,
        L_dimensionless=L,strip_half_width=interval(rho),strip_count=4,terms=terms,
        input_file=source.relative_to(root).as_posix(),input_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        output_norm_squared=interval(nw),gaussian_overlap_tail=interval(tail),
        strip_overlap_bound=interval(strips),product_overlap_bound=interval(product),
        total_omission_bound=interval(error),
        hashes={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in
                ['code/residual_gaussian_overlap_bounded.py','code/residual_overlap_rectangle.py','code/chirp_residual_certificate.py','code/gaussian_outer_bound.py']},
        versions=dict(python=platform.python_version(),python_flint=flint.__version__))
    (base/'gaussian-bounded-input.json').write_text(json.dumps(metadata,indent=2))
    total=acb(0);finite=True;start=time.perf_counter()
    with log.open('x') as f:
        for i,(ya,yb) in enumerate(intervals):
            calls=dict(outer=0,inner=0,rejected_outer=0);t=time.perf_counter()
            def outer(y,analytic):
                calls['outer']+=1
                if not outer_guard(y):
                    calls['rejected_outer']+=1
                    return acb('nan')
                if analytic or y.real.rad()>arb('1e-25') or y.imag.rad()>arb('1e-25'):
                    return gaussian_outer_bound(y,n,normw,terms)
                pw=acb_poly([p(y) for p in pb])
                def inner(x,analytic_x):
                    calls['inner']+=1
                    return pw(x)*(-x*x-y*y).exp()*phase_product(scale*x,scale*y,terms)/pi
                return sum((acb.integral(inner,j,j+3,abs_tol=arb('1e-10'),rel_tol=arb('1e-10'),
                       eval_limit=3000,depth_limit=40) for j in [-6,-3,0,3]),acb(0))
            value=acb.integral(outer,ya,yb,abs_tol=arb('1e-8'),rel_tol=arb('1e-8'),eval_limit=3000,depth_limit=40)
            row=dict(segment=i,lower=interval(ya),upper=interval(yb),finite=bool(value.is_finite()),
                     seconds=time.perf_counter()-t,calls=calls)
            if value.is_finite():
                row.update(real=interval(value.real),imag=interval(value.imag));total+=value
            else:finite=False
            f.write(json.dumps(row)+'\n');f.flush()
            print(json.dumps(dict(segment=i,finite=row['finite'],seconds=row['seconds'],calls=calls,value=str(value))),flush=True)
    metadata.update(status='full Gaussian overlap enclosed' if finite else 'nonfinite segment; full overlap unresolved',
        seconds=time.perf_counter()-start,segments=len(intervals),
        scope='One of five F overlaps only. No full word certificate.')
    if finite:
        lower=total.real-error
        metadata.update(compact_real=interval(total.real),compact_imag=interval(total.imag),
            full_real_lower_bound=interval(lower),deficit_upper_bound=interval(nw-lower),
            target_eta='0.0001',target_proved=bool(lower>nw-arb('0.0001')))
    out.write_text(json.dumps(metadata,indent=2))
    print(json.dumps(dict(status=metadata['status'],target_proved=metadata.get('target_proved'))),flush=True)

if __name__=='__main__':main()
