"""Full-plane F-step overlap with explicitly bounded strips and tails.

Results are one-shot. Every integration segment is saved before the next.
The first run targets only right step 1 of saved case 06, not the full word.
"""
from pathlib import Path
import sys,json,hashlib,time,platform
root=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(root/'tmp/arb-runtime'))
import numpy as np
import flint
from flint import arb,acb,ctx
from chirp_residual_certificate import complex_exact,norm_squared,interval
from residual_overlap_rectangle import integrate_rectangle

def coordinate_moment(c,power,transpose=False):
    if transpose:c=c.T
    n=c.shape[0]
    a=[[complex_exact(z) for z in row] for row in c]
    for _ in range(power):
        old=len(a);b=[]
        for i in range(old+1):
            lo=(arb(i)/2).sqrt() if i else arb(0)
            hi=(arb(i+1)/2).sqrt()
            b.append([(lo*a[i-1][j] if i else 0)+(hi*a[i+1][j] if i+1<old else 0)
                      for j in range(n)])
        a=b
    return sum((z.real*z.real+z.imag*z.imag for row in a for z in row),arb(0))

def square_tail(c,L,power):
    return (coordinate_moment(c,power)+coordinate_moment(c,power,True))/arb(L)**(2*power)

def main():
    base=root/'results/stage4/residual-certificate'
    out=base/'full-right-1.json';log=base/'full-right-1-segments.jsonl'
    if out.exists() or log.exists():raise RuntimeError('Do not overwrite or restart a running study')
    ctx.prec=768;L=20;rho=arb('1e-8');terms=10;power=32
    source=root/'results/stage4/flip-domain/coefficients-06.npz'
    with np.load(source,allow_pickle=False) as data:v,w=data['right'][0:2]
    nv,nw=norm_squared(v),norm_squared(w);norm_product=(nv*nw).sqrt()
    start=time.perf_counter();tv,tw=square_tail(v,L,power),square_tail(w,L,power)
    tails=(tv*tw).sqrt()
    r=arb(2*v.shape[0]-1)/(2*v.shape[0])
    mehler=r**(-(v.shape[0]-1))/(arb.pi()*(1-r*r)).sqrt()
    poles=[arb(k)*arb.pi().sqrt() for k in range(-19,20,2)
           if abs(arb(k)*arb.pi().sqrt()) < L]
    strips=2*rho*len(poles)*mehler*norm_product
    first=(arb(L)*arb.pi().sqrt()-(2*terms+1)*arb.pi()).exp()
    products=4*first/((1-(-2*arb.pi()).exp())*(1-first))*norm_product
    error=tails+strips+products
    metadata=dict(status='running',index=6,n=128,branch='right',step=1,L_dimensionless=L,
        strip_half_width=interval(rho),strip_count=len(poles),moment_power=power,terms=terms,
        precision_bits=ctx.prec,input_file=source.relative_to(root).as_posix(),
        input_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        input_norm_squared=interval(nv),output_norm_squared=interval(nw),
        input_square_tail_bound=interval(tv),output_square_tail_bound=interval(tw),
        overlap_tail_bound=interval(tails),strip_overlap_bound=interval(strips),
        product_overlap_bound=interval(products),omission_bound=interval(error),
        versions=dict(python=platform.python_version(),python_flint=flint.__version__),
        hashes={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in
                ['code/residual_full_overlap.py','code/residual_overlap_rectangle.py','code/chirp_residual_certificate.py']})
    (base/'full-right-1-input.json').write_text(json.dumps(metadata,indent=2))
    print(json.dumps(dict(event='bounds',omission=error.str(15),tail=tails.str(15),strips=strips.str(15))),flush=True)
    intervals=[];a=arb(-L)
    for p in poles:
        intervals.append((a,p-rho));a=p+rho
    intervals.append((a,arb(L)))
    total=acb(0);finite=True
    with log.open('x') as f:
        for i,(a,b) in enumerate(intervals):
            t=time.perf_counter()
            val,calls=integrate_rectangle(v,w,-L,L,a,b,terms=terms)
            row=dict(segment=i,lower=interval(a),upper=interval(b),finite=bool(val.is_finite()),
                     seconds=time.perf_counter()-t,calls=calls)
            if val.is_finite():
                row.update(real=interval(val.real),imag=interval(val.imag));total+=val
            else:finite=False
            f.write(json.dumps(row)+'\n');f.flush()
            print(json.dumps(dict(segment=i,finite=row['finite'],seconds=row['seconds'],
                 value=str(val),calls=calls)),flush=True)
    metadata.update(status='full overlap enclosed' if finite else 'nonfinite segment; no full enclosure',
        segments=len(intervals),seconds=time.perf_counter()-start,
        scientific_status='Only one of five F-step overlaps targeted; no continuum pentagon witness.')
    if finite:
        lower=total.real-error;upper=total.real+error
        metadata.update(compact_real=interval(total.real),compact_imag=interval(total.imag),
            full_real_lower_bound=interval(lower),full_real_upper_bound=interval(upper),
            deficit_upper_bound=interval(nw-lower),target_eta='0.0001',
            target_proved=bool(lower>nw-arb('0.0001')))
    out.write_text(json.dumps(metadata,indent=2))
    print(json.dumps(dict(status=metadata['status'],target_proved=metadata.get('target_proved'))),flush=True)

if __name__=='__main__':main()
