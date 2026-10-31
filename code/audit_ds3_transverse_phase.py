"""Independent saved-arithmetic and symbolic transverse checks; no interval claim."""
import hashlib,json
from pathlib import Path
from decimal import Decimal as D,localcontext
from fractions import Fraction as F
from math import comb
import sympy as sp
from ds3_peak_normalization import gamma_terms
from ds3_full_tail_bound import parameters
from audit_ds3_peak_limit_refined import machin


def symbolic():
    M,u,v,p,q,r,L=sp.symbols('M u v p q r L',real=True)
    a=M+u+v;b=M+u;c=M;d=M+v+r
    hs=[(1+sp.I*x)/2 for x in [a,q,c,p]];h=(1+sp.I*b)/2;hp=(1+sp.I*d)/2
    _,X=parameters(h,hp,hs);_,Y=parameters(1-h,hp,hs)
    same_top=all(sp.simplify(x-y)==0 for x,y in zip(X[0][0],Y[0][0]))
    same_bottom=sorted(map(str,X[0][1]))==sorted(map(str,Y[0][1]))
    terms=gamma_terms(b,d,a,q,c,p,'peak_phase');phase=0;small={};large={}
    for arg,w in terms.items():
        slope=sp.expand(sp.im(arg)).coeff(M)
        if slope==0:small[arg]=w;continue
        large[arg]=w;shift=sp.expand(arg-sp.I*slope*M)
        expr=(sp.I*slope*M+shift-sp.Rational(1,2))*(L+sp.log(abs(slope))+sp.I*sp.sign(slope)*sp.pi/2)-sp.I*slope*M
        expr+=sp.bernoulli(2,shift)/(2*sp.I*slope*M)
        phase+=w*sp.im(sp.expand(expr))
    phase=sp.factor(phase);expected=-r*L/2-r*(3*r+4*v)/(8*M)
    A=[(1+sp.I*x)/2 for x in [v+q,p-u,-u-p,v-q]];z=-sp.I*r/2
    # Remaining gamma terms outside the declared D_infinity phase must be conjugate pairs.
    declared={}
    for arg,w in [(sp.conjugate(A[0]),1),(sp.conjugate(A[1]),1),
                  (sp.conjugate(A[0])+z,-1),(sp.conjugate(A[1])+z,-1),(1-z,-1)]:
        arg=sp.expand(arg);declared[arg]=declared.get(arg,0)+w
    remainder=dict(small)
    for arg,w in declared.items():remainder[arg]=remainder.get(arg,0)-w
    pairs=all(w==remainder.get(sp.expand(sp.conjugate(arg)),0) for arg,w in remainder.items())
    return dict(shared_top=same_top,shared_bottom=same_bottom,large_phase_correct=sp.simplify(phase-expected)==0,
                small_phase_correct=pairs,large_phase=str(phase),
                small_terms=[dict(argument=str(k),weight=w) for k,w in small.items()],
                large_terms=[dict(argument=str(k),weight=w) for k,w in large.items()])


def audit(root):
    data=json.loads((root/'results/stage4/ds3-transverse-phase/study.json').read_text());errors=[];count=0
    def req(ok,msg):
        nonlocal count
        count+=1
        if not ok:errors.append(msg)
    def z(x):return tuple(D(v) for v in x)
    def mul(a,b):return a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]
    def norm(a):return (a[0]*a[0]+a[1]*a[1]).sqrt()
    def close(a,b,msg):req(abs(a-b)<D('1e-75')*max(D(1),abs(b)),msg)
    with localcontext() as ctx:
        ctx.prec=130;pi,_,_=machin()
        req(data['status']=='passed' and len(data['rows'])==16,'16 saved rows')
        for i,row in enumerate(data['rows']):
            D0=z(row['D0']);E0=z(row['E0']);S=z(row['S']);val=z(row['value'])
            finite=2*mul(D0,S)[0];close(finite,D(row['finite_leading']),f'{i}: finite leading')
            err=norm((val[0]-finite,val[1]));close(err,D(row['finite_error']),f'{i}: error')
            RF,RY,BG=map(D,[row['RF'],row['RY'],row['BG']])
            bound=2*norm(D0)*(RY+RF*(norm(S)+RY))+2*norm(E0)*BG*(norm(S)+RY)
            close(bound,D(row['finite_bound']),f'{i}: full bound')
            req(err<=bound+D(row['tail']) and row['finite_passed'],f'{i}: finite inequality')
            b=D(row['M'])+D(row['u']);d0=D(row['M'])+D(row['v']);r=D(row['r'])
            ma=min(d0+r,b-r/2);mb=min(d0+r/2,b)
            ya,yb=map(z,row['coefficients'][2:])
            close((pi*pi/6-1)/(1+ma*ma),RF,f'{i}: RF')
            close((pi*pi/6-1)*(norm(ya)/(1+ma*ma)+norm(yb)/(1+mb*mb)),RY,f'{i}: RY')
            poly=list(map(F,row['details']['omega'][0]['details']['groups'][1]['certificate']['coefficients']))
            K=max(abs(x) for x in poly[1:]);N=2+(K.numerator+K.denominator-1)//K.denominator
            req(N==row['Nstar'],f'{i}: auxiliary threshold')
            ascending=list(reversed(poly))
            req(all(sum(ascending[j]*comb(j,k)*N**(j-k) for j in range(k,11))>0 for k in range(11)),f'{i}: exact translated positivity')
            close((3*N+3)*D(row['Bb']),BG,f'{i}: companion bound')
            req(D(row['scattering_identity_error'])<D('1e-75'),f'{i}: scattering coefficient')
            req(abs(D(row['nonoscillatory_real']))<D('1e-75'),f'{i}: nonoscillatory cancellation')
            close(norm(z(row['T'])),1/(2*pi*abs(r)),f'{i}: exact scattering modulus')
        for path,digest in data['hashes'].items():req(hashlib.sha256((root/path).read_bytes()).hexdigest()==digest,path)
    sym=symbolic()
    for k in ['shared_top','shared_bottom','large_phase_correct','small_phase_correct']:req(sym[k],k)
    return dict(status='failed' if errors else 'passed',checks=count,errors=errors,symbolic=sym,
                scope='Independent Decimal bound/residual replay, exact rational threshold and symbolic Gamma phase; not a uniform transverse or interval arithmetic certification.')


if __name__=='__main__':
    root=Path(__file__).resolve().parents[1];dest=root/'results/stage4/ds3-transverse-phase/audit.json'
    if dest.exists():raise RuntimeError('Preserve prior audit')
    result=audit(root);result['audit_source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print(result['status'],result['checks'],result['errors'])
    raise SystemExit(result['status']!='passed')
