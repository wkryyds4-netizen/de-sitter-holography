"""Independent finite-bound and polynomial-threshold replay for the collision peak."""
import hashlib,json
from pathlib import Path
from decimal import Decimal as D,localcontext
from fractions import Fraction as F
from math import comb
import sympy as sp
from audit_ds3_crossing_dictionary import PI


def audit(root):
    data=json.loads((root/'results/stage4/ds3-peak-limit/study.json').read_text());errors=[];checks=0
    def req(ok,msg):
        nonlocal checks
        checks+=1
        if not ok:errors.append(msg)
    def norm(z):return sum(D(v)**2 for v in z).sqrt()
    def relative(a,b):return abs(a-b)<D('1e-85')*max(abs(a),abs(b),D('1e-500'))
    symbolic=data['symbolic'];terms={sp.sympify(r['argument'],locals={k:sp.Symbol(k,real=True) for k in ['b','d','a','q','c','p']}):r['weight'] for r in symbolic['coefficient_terms']}
    req(all(w==terms.get(sp.expand(sp.conjugate(z)),0) for z,w in terms.items()),'full coefficient Gamma conjugate pairs')
    req(all(symbolic[k] for k in ['coefficient_conjugate_pairs','small_phase_cancels','phase_expansion_correct']),'symbolic checks')
    req(len(data['rows'])==12 and data['status']=='passed','twelve full-kernel samples')
    with localcontext() as ctx:
        ctx.prec=150
        for row in data['rows']:
            M,u,v,p,q=map(D,[row['M'],row['u'],row['v'],row['p'],row['q']]);b=M+u;d=M+v;A=D(row['Amax'])
            coeff=list(map(F,row['details']['omega'][0]['details']['groups'][1]['certificate']['coefficients']))
            mx=max(map(abs,coeff[1:]));N=2+(mx.numerator+mx.denominator-1)//mx.denominator
            req(N==row['Nstar'],'explicit Cauchy threshold')
            poly=coeff[::-1]
            shifted=[sum(poly[j]*comb(j,k)*N**(j-k) for j in range(k,len(poly))) for k in range(len(poly))]
            req(all(x>0 for x in shifted),'all derivative coefficients positive at threshold')
            B=(((PI*b).exp()-(-PI*b).exp())/(2*PI*b)).sqrt()*(PI*(2*b*A+A*A)/b).exp()
            big=(3*N+3)*B;req(relative(big,D(row['large_series_bound'])),'uniform large-series bound')
            E=norm(row['suppressed_coefficient']);req(relative(E*big,D(row['suppressed_term_bound'])),'suppressed coefficient times whole-series bound')
            RF=(PI*PI/6-1)/(1+min(b,d)**2);req(relative(RF,D(row['RF'])),'small-series error')
            D0=[D(x) for x in row['leading_coefficient']];phase_distance=norm([4*PI*D0[0]-1,4*PI*D0[1]])
            RY=D(row['RY']);q0=row['q0']
            bound=RY/(2*PI)+((RF+phase_distance)/(2*PI)+2*E*big)*(norm(q0)+RY)
            req(relative(bound,D(row['leading_bound'])),'finite-M whole-kernel bound')
            leading=D(q0[0])/(2*PI);req(relative(leading,D(row['leading'])),'leading real digamma value')
            error=norm([D(row['value'][0])-leading,D(row['value'][1])])
            req(abs(error-D(row['leading_error']))<D('1e-90'),'complex leading residual')
            req(error<=bound+D(row['numerical_tail']) and row['passed'],'observed residual and declared tail')
            theta2=(-p*p*v+q*q*u+u*u*v-u*v*v+u-v)/4
            req(abs(theta2-D(row['phase_coefficient']))<D('1e-90'),'phase coefficient')
            req(relative(D(row['phase'])*M*M,D(row['phase_scaled'])),'scaled phase')
            req(D(row['unit_modulus_residual'])<D('1e-90') and D(row['coefficient_conjugacy_residual'])<D('1e-90'),'exact coefficient identities sampled')
            approx=M.ln()/PI+D(row['constant'])+(u+v)/(2*PI*M)
            req(relative(approx,D(row['expansion'])),'logarithm constant and inverse mass term')
            err=norm([D(row['value'][0])-approx,D(row['value'][1])])
            req(abs(err-D(row['expansion_error']))<D('1e-90'),'expansion residual')
            req(relative(err*M*M/M.ln(),D(row['scaled_expansion_error'])),'scaled remainder')
        req(D(data['rows'][0]['leading_bound'])>D('1e4'),'retain loose low-mass bound')
    for path,sha in data['hashes'].items():req(hashlib.sha256((root/path).read_bytes()).hexdigest()==sha,'hash '+path)
    return dict(status='passed' if not errors else 'failed',checks=checks,errors=errors,
                scope='Saved arithmetic, exact rational polynomial thresholds and Gamma-phase certificate. No interval rounding or mass-integral convergence certification.')


if __name__=='__main__':
    root=Path(__file__).resolve().parents[1];dest=root/'results/stage4/ds3-peak-limit/audit.json'
    if dest.exists():raise RuntimeError('Preserve prior audit')
    data=audit(root);data['audit_source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    dest.write_text(json.dumps(data,indent=2),encoding='utf-8');print(json.dumps(data,indent=2));raise SystemExit(data['status']!='passed')
