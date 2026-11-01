"""Ball-arithmetic Gaussian moments for saved finite Hermite vectors.

No quadrature is used for chirp steps. Bounds on the other five overlaps
remain hypotheses; the conditional margin is not a continuum witness.
Install python-flint==0.9.0, or use the project-local tmp/arb-runtime.
"""
from pathlib import Path
import sys, json, hashlib, time, platform
root = Path(__file__).resolve().parents[1]
if (root/'tmp/arb-runtime').is_dir():
    sys.path.insert(0, str(root/'tmp/arb-runtime'))
import numpy as np
import flint
from flint import arb, acb, acb_mat, acb_poly, ctx, fmpq

def real_exact(x):
    a,b = float(x).as_integer_ratio()
    return arb(fmpq(a,b))

def complex_exact(z):
    return acb(real_exact(z.real), real_exact(z.imag))

def norm_squared(c):
    return sum((real_exact(z.real)**2+real_exact(z.imag)**2 for z in c.flat), arb(0))

def norm_difference(a,b):
    # Subtract after exact lifting, not in numpy binary64.
    return sum(((real_exact(x.real)-real_exact(y.real))**2+
                (real_exact(x.imag)-real_exact(y.imag))**2
                for x,y in zip(a.flat,b.flat)), arb(0)).sqrt()

def polynomial_matrix(c):
    n = c.shape[0]
    ps = [acb_poly([1])]
    if n > 1:
        ps.append(acb_poly([0,arb(2).sqrt()]))
    for k in range(1,n-1):
        ps.append(acb_poly([0,(arb(2)/(k+1)).sqrt()])*ps[-1]
                  -(arb(k)/(k+1)).sqrt()*ps[-2])
    t = acb_mat(n,n,[ps[j][i] for i in range(n) for j in range(n)])
    a = acb_mat(n,n,[complex_exact(z) for z in c.flat])
    return t*a*t.transpose()

def chirp_overlap(v,w,sign):
    """<w, exp(i sign z t) v> in normalized oscillator coordinates."""
    assert sign in [-1,1] and v.shape == w.shape
    n = v.shape[0]
    a,b = polynomial_matrix(v),polynomial_matrix(w)
    p = [acb_poly([a[i,j] for j in range(n)]) for i in range(n)]
    q = [acb_poly([b[i,j].conjugate() for j in range(n)]) for i in range(n)]
    prod = [acb_poly([]) for _ in range(2*n-1)]
    for i in range(n):
        for j in range(n):
            prod[i+j] += p[i]*q[j]
    size = 2*n-1
    diagonal = arb(2)/5
    off_diagonal = acb(0,arb(sign)/5)
    mom = [[acb(0) for _ in range(size)] for _ in range(size)]
    mom[0][0] = acb(1)
    for j in range(2,size,2):
        mom[0][j] = (j-1)*diagonal*mom[0][j-2]
    val = acb(0)
    for i in range(size):
        for j in range(size):
            if i:
                mom[i][j] = ((i-1)*diagonal*mom[i-2][j] if i>1 else 0) + (
                    j*off_diagonal*mom[i-1][j-1] if j else 0)
            val += prod[i][j]*mom[i][j]
    return 2/arb(5).sqrt()*val

def interval(x):
    def pair(q):
        q=q.fmpq()
        return [str(q.p),str(q.q)]
    return dict(lower=pair(x.lower()),upper=pair(x.upper()),display=x.str(22))

def replay_case(index,precision):
    ctx.prec = precision
    path = root/f'results/stage4/flip-domain/coefficients-{index:02d}.npz'
    with np.load(path,allow_pickle=False) as data:
        left,right = data['left'],data['right']
    records=[];chirp_errors=[]
    for step,sign in [(2,1),(4,-1)]:
        start=time.perf_counter();v,w=right[step-1],right[step]
        ov=chirp_overlap(v,w,sign)
        e2=norm_squared(v)+norm_squared(w)-2*ov.real
        assert e2 > 0 and e2.rad() < arb('1e-24'), str(e2)
        e=e2.sqrt();chirp_errors.append(e)
        records.append(dict(step=step,sign=sign,overlap_real=interval(ov.real),
            overlap_imag=interval(ov.imag),residual_squared=interval(e2),
            residual=interval(e),seconds=time.perf_counter()-start))
    residual=norm_difference(left[-1],right[-1])
    hypotheses=[]
    for branch,vec,steps in [('left',left,[1,2]),('right',right,[1,3,5])]:
        for step in steps:
            nv,nw=norm_squared(vec[step-1]),norm_squared(vec[step])
            hypotheses.append(dict(branch=branch,step=step,
                input_norm_squared=interval(nv),output_norm_squared=interval(nw)))
    conditional=[]
    for eta_text in ['0.000001','0.00001','0.0001']:
        eta=arb(eta_text);errors=[]
        for vec,steps in [(left,[1,2]),(right,[1,3,5])]:
            for step in steps:
                e2=norm_squared(vec[step-1])-norm_squared(vec[step])+2*eta
                assert e2 > 0
                errors.append(e2.sqrt())
        margin=residual-sum(chirp_errors,arb(0))-sum(errors,arb(0))
        conditional.append(dict(eta=eta_text,margin=interval(margin),positive=bool(margin>0),
            hypothesis='Each of the five uncomputed exact F-step overlaps has real part >= output norm squared - eta.'))
    return dict(index=index,n=left.shape[1],precision_bits=precision,input_file=path.relative_to(root).as_posix(),
        input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),chirp_steps=records,
        final_vector_separation=interval(residual),unverified_overlap_targets=hypotheses,
        conditional_margins=conditional)

def main():
    out=root/'results/stage4/residual-certificate/chirp.json'
    if out.exists():raise RuntimeError('Do not overwrite')
    out.parent.mkdir(parents=True,exist_ok=True)
    ctx.prec=512
    f=np.ones((1,1),complex)
    controls=[]
    for sign in [-1,1]:
        val=chirp_overlap(f,f,sign)
        expected=2/arb(5).sqrt()
        assert val.real.overlaps(expected) and val.imag.contains(0)
        # Independent one-dimensional analytic Gaussian integral, including t integration.
        direct=acb.integral(lambda z,_:(-arb(5)*z*z/4).exp(),-12,12,
                            abs_tol=arb('1e-40'),rel_tol=arb('1e-40'))/arb.pi().sqrt()
        tail=4*(-arb(5)*144/4).exp()/(arb(5)*12*arb.pi().sqrt())
        assert abs(direct.real-expected) < tail+arb('1e-38')
        controls.append(dict(sign=sign,moment_real=interval(val.real),
            independent_integral_real=interval(direct.real),tail_bound=interval(tail)))
    rows=[]
    for index,precision in [(3,512),(6,768),(6,1024)]:
        row=replay_case(index,precision);rows.append(row)
        print(json.dumps(dict(index=index,precision=precision,
            chirp=[x['residual']['display'] for x in row['chirp_steps']],
            conditional=[(x['eta'],x['margin']['display']) for x in row['conditional_margins']])),flush=True)
    inputs=['code/chirp_residual_certificate.py','results/stage4/flip-domain/evidence.json']
    out.write_text(json.dumps(dict(status='chirp steps certified; F-step overlaps uncomputed',
        controls=controls,rows=rows,
        versions=dict(python=platform.python_version(),python_flint=flint.__version__,numpy=np.__version__),
        hashes={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in inputs},
        scope='Exact lifting of saved binary64 coefficients. Gaussian moments enclosed by Arb; no continuum pentagon lower bound without the five overlap hypotheses.'),indent=2))

if __name__=='__main__':main()
