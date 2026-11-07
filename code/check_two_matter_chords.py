"""Channel decomposition, independent Wick check, and relative-mass checks."""
from pathlib import Path
from itertools import combinations
import json,hashlib,platform
import sympy as s
import mpmath as mp
from two_matter_chords import TwoMatter,basis
from singlet_coefficient import pairings
from singlet_multimark import crossing
ROOT=Path(__file__).resolve().parent.parent;OUT=ROOT/'results/susskind-two-matter'

def symbolic():
    r,chi=s.symbols('r chi',positive=True)
    G=s.Matrix([[1,r,r*r],[r,1,r],[r*r,r,1]])
    X=chi*s.Matrix([[1,r,r*r],[r,r*r,r],[r*r,r,1]])
    primary=s.Matrix([-r,1+r*r,-r])/s.sqrt(1-r**4)
    outer=s.Matrix([1,0,1]);v=s.Matrix([1,-2,1])
    a=-2*(1-r*r)/s.sqrt(1-r**4);c=(1-r)**2/(1+r*r)
    tests={}
    def check(name,expr):
        ok=all(s.simplify(x)==0 for x in expr) if isinstance(expr,s.MatrixBase) else s.simplify(expr)==0
        assert ok,name;tests[name]=True
    check('primary_direct_norm',(primary.T*G*primary)[0]-1)
    check('primary_crossed_overlap',(primary.T*X*primary)[0]+chi*r*r)
    check('primary_external_orthogonality',(primary.T*G*outer)[0])
    check('primary_annihilation',s.Matrix([[1,r,r*r],[r*r,r,1]])*primary)
    check('v_decomposition',v-c*outer-a*primary)
    two=s.factor((v.T*(G+X)*v)[0].subs(chi,r*r))
    parts=[4*(1-r)**2,2*(1-r)**4,4*(1-r*r)**2]
    check('two_sector_sum',two-parts[1]-parts[2])
    total=s.factor(sum(parts));expected=2*(1-r)**2*(5+2*r+3*r*r)
    check('total_exact_coefficient',total-expected)
    # Independent cyclic Wick enumeration of F^2; every pairing retained.
    q=s.symbols('q',positive=True)
    def moment(word):
        hs=tuple(i for i,z in enumerate(word) if z=='H');xs=tuple(i for i,z in enumerate(word) if z=='X')
        out=0
        for hp in pairings(hs):
            for xp in pairings(xs):
                e=sum(crossing(a,b) for a,b in combinations(hp,2))
                e+=2*sum(crossing(a,b) for a in hp for b in xp)
                e+=4*sum(crossing(a,b) for a,b in combinations(xp,2))
                out+=q**e
        return s.expand(out)
    words=[('HXX',1),('XHX',-2),('XXH',1)]
    wick=s.expand(sum(a*b*moment(u+w) for u,a in words for w,b in words))
    check('independent_six_letter_Wick',wick-total.subs(r,q*q))
    h=s.symbols('h',positive=True)
    slopes=[s.limit(4*z.subs(r,s.exp(-2*h))/h**2,h,0) for z in parts]
    assert slopes==[64,0,256]
    # Exact d1=d2=2 Jacobi basis and free mass matrix.
    x=s.symbols('x',positive=True);weight=x**3*(1-x)**3
    pol=[s.jacobi(k,3,3,2*x-1) for k in range(5)]
    norms=[s.integrate(s.expand(weight*p*p),(x,0,1)) for p in pol]
    mass=[]
    for j in range(5):
        row=[]
        for k in range(5):
            norm=s.integrate(s.expand(weight*pol[j]*pol[k]),(x,0,1))/s.sqrt(norms[j]*norms[k])
            check(f'jacobi_orthogonality_{j}_{k}',norm-int(j==k))
            value=s.integrate(s.expand(x**2*(1-x)**2*pol[j]*pol[k]),(x,0,1))/s.sqrt(norms[j]*norms[k])
            row.append(str(s.simplify(value)))
        mass.append(row)
    # Jacobi differential eigenproblem gives relative Casimir, not mass.
    for k,p in enumerate(pol):
        L=-x*(1-x)*s.diff(p,x,2)-(4-8*x)*s.diff(p,x)
        check(f'relative_Casimir_{k}',L-k*(k+7)*p)
    assert s.sympify(mass[0][2])!=0
    return dict(checks=tests,wick_polynomial=str(wick),exact_total=str(total),channel_norm_numerators=[str(z) for z in parts],C_t2_over_Js2h_limits=[str(z) for z in slopes],free_mass_matrix=mass)

def numerical(h,d1,d2,prec):
    mp.mp.dps=prec;A=TwoMatter(h,d1,d2);B=TwoMatter(h,d2,d1);rows=[]
    for k in range(5):
        u=A.primary(k);v=B.primary(k);G=A.matrix(k);X=A.matrix(k,True)
        overlap=(u.T*X*v)[0]
        pred=(-1)**k*A.q**(mp.mpf(k*(k-1))/2+k*(A.d1+A.d2)+A.d1*A.d2)
        err=abs(overlap-pred);assert err<mp.mpf('1e-60')
        mid=basis(k).index((0,k,0));psi=(u.T*G)[mid]
        expected=mp.fprod(A.qi(j) for j in range(1,k+1))
        expected*=mp.fprod((1-A.q**(2*A.d1+j))*(1-A.q**(2*A.d2+j))/(1-A.q**(2*A.d1+2*A.d2+k-1+j)) for j in range(k))
        assert abs(psi*psi-expected)<mp.mpf('1e-60')
        # k=1, identical matter: direct+cross is the physical two-open-X form.
        if k==1 and d1==d2:
            rr=A.q**A.d1
            vcomm=mp.matrix([1,-2,1]) # basis(1)=(R,middle,L)
            got=(vcomm.T*(G+X)*vcomm)[0]
            target=2*(1-rr)**2*(3+2*rr+3*rr*rr)
            assert abs(got-target)<mp.mpf('1e-60')
        rows.append(dict(k=k,overlap=str(overlap),formula=str(pred),overlap_abs_error=str(err),middle_projection_squared=str(psi*psi)))
    return dict(h=h,d1=d1,d2=d2,dps=prec,rows=rows)

def main():
    OUT.mkdir(parents=True,exist_ok=True);target=OUT/'validation.json'
    if target.exists():raise FileExistsError(target)
    exact=symbolic();cases=[];ref=[]
    for h,d1,d2 in [('0.05','2','2'),('0.01','2','2'),('0.02','1.3','2.1')]:
        a=numerical(h,d1,d2,85);b=numerical(h,d1,d2,115)
        worst=max(abs(mp.mpf(x[key])/mp.mpf(y[key])-1) for x,y in zip(a['rows'],b['rows']) for key in ['overlap','middle_projection_squared'])
        assert worst<mp.mpf('1e-55')
        cases.append(b);ref.append(dict(h=h,d1=d1,d2=d2,max_relative_difference=str(worst)))
        print(h,d1,d2,'passed',flush=True)
    paths=[Path(__file__),ROOT/'code/two_matter_chords.py',ROOT/'code/singlet_coefficient.py',ROOT/'code/singlet_multimark.py']
    data=dict(status='passed',versions=dict(python=platform.python_version(),mpmath=mp.__version__,sympy=s.__version__),symbolic=exact,cases=cases,precision_refinement=ref,
        scope='Exact first-grade channel proof plus independent full Wick sum. Numerical k<=4 checks do not prove the general conjectured LS overlap formula.',
        input_sha256={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    target.write_text(json.dumps(data,indent=2)+'\n')
    print('mass 0,2:',exact['free_mass_matrix'][0][2]);print('saved',target)
if __name__=='__main__':main()
