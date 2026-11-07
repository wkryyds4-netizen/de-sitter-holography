"""Casimir form of the leading ordering-odd kernel and translation compatibility."""
from pathlib import Path
import sympy as s
import json,hashlib,platform
ROOT=Path(__file__).resolve().parent.parent;OUT=ROOT/'results/susskind-ordering-metric'

def main():
    OUT.mkdir(parents=True,exist_ok=True);target=OUT/'exact.json'
    if target.exists():raise FileExistsError(target)
    tests={}
    def check(name,res):
        assert s.simplify(res)==0,name;tests[name]=True
    d,e,a,ap,b,bp=s.symbols('d e a ap b bp',positive=True)
    def elements(d,a,ap):return [d*(1+a*ap)/(a+ap),s.I*d*(ap-a)/(a+ap),d*(1-a*ap)/(a+ap)]
    E1,B1,P1=elements(d,a,ap);E2,B2,P2=elements(e,b,bp)
    kernel=2*(E1*E2-B1*B2-P1*P2)
    expected=2*d*e*(2*a*ap+2*b*bp+(ap-a)*(bp-b))/((a+ap)*(b+bp))
    check('coherent_Casimir_equals_commutator_covariance',kernel-expected)
    x=s.Symbol('x',positive=True);f=s.Function('f')(x);M=s.Function('M')(x)
    def L(v):return -x*(1-x)*s.diff(v,x,2)-(2*d-2*(d+e)*x)*s.diff(v,x)
    A=lambda v:L(v)+2*d*e*v
    check('general_multiplication_commutator',A(M*f)-M*A(f)+2*x*(1-x)*s.diff(M,x)*s.diff(f,x)
          +(x*(1-x)*s.diff(M,x,2)+(2*d-2*(d+e)*x)*s.diff(M,x))*f)
    # Conjugated relative Casimir from the previously verified tensor coordinate map.
    w0=x**(d-s.Rational(1,2))*(1-x)**(e-s.Rational(1,2))
    T=lambda v:-s.diff(x*(1-x)*s.diff(v,x),x)+((d-s.Rational(1,2))**2/x+(e-s.Rational(1,2))**2/(1-x)+s.Rational(1,4))*v
    check('relative_Casimir_conjugation',s.expand(T(w0*f)/w0-s.Rational(1,4)*f-d*(d-1)*f-e*(e-1)*f-A(f)))
    weight=504*x**3*(1-x)**5 # beta(4,6) probability density
    check('beta_normalization',s.integrate(weight,(x,0,1))-1)
    pol=[];norms=[];ak=[]
    for k in range(6):
        pk=s.jacobi(k,5,3,2*x-1);nk=s.integrate(s.expand(weight*pk**2),(x,0,1))
        pol.append(pk/s.sqrt(nk));norms.append(nk);ak.append(k*(k+9)+12)
        check(f'Jacobi_eigen_{k}',s.expand(A(pk).subs({d:2,e:3})-ak[-1]*pk))
    mat=s.zeros(6);mass=1/(x*(1-x))
    for j in range(6):
      for k in range(6):
        check(f'orthogonality_{j}_{k}',s.integrate(s.expand(weight*pol[j]*pol[k]),(x,0,1))-int(j==k))
        mat[j,k]=s.simplify(s.integrate(s.expand(504*x**2*(1-x)**4*pol[j]*pol[k]),(x,0,1)))
    check('mass_00',mat[0,0]-s.Rational(24,5));check('mass_01',mat[0,1]+s.sqrt(66)/10)
    Am=s.diag(*ak);root=s.diag(*[s.sqrt(z) for z in ak]);naive=root*mat*root.inv()
    check('mass_Casimir_defect_01',(Am*mat-mat*Am)[0,1]-s.sqrt(66))
    assert s.simplify(naive[0,1]-naive[1,0])!=0
    # Exact norm under the undressed multiplication phase, in the A-weighted space.
    state=(pol[0]+s.I*pol[1])/s.sqrt(2);fstar=s.conjugate(state);p=x*(1-x)
    linear=s.I*p*s.diff(mass,x)*(fstar*s.diff(state,x)-s.diff(fstar,x)*state)
    quad=p*s.diff(mass,x)**2*fstar*state
    c1=s.simplify(s.integrate(s.cancel(weight*linear),(x,0,1)))
    c2=s.simplify(s.integrate(s.cancel(weight*quad),(x,0,1)))
    check('phase_norm_linear',c1-s.sqrt(66));norm=s.Rational(ak[0]+ak[1],2)
    rows=[[str(mat[j,k]) for k in range(6)] for j in range(6)]
    payload=dict(status='passed',checks=tests,relative_labels=ak,mass_matrix=rows,
        naive_mass_in_A_orthonormal_basis=[[str(s.simplify(naive[j,k])) for k in range(2)] for j in range(2)],
        undressed_phase_norm=dict(state='(p0+i p1)/sqrt(2), p0=1, p1=(10x-4)/sqrt(24/11)',constant=str(norm),linear=str(c1),quadratic=str(c2)),
        ward_defect_01=str(s.simplify((root*mat-mat*root)[0,1])),
        general_operator='A = -x(1-x) d_x^2 - [2d-2(d+e)x] d_x + 2de',
        environment=dict(python=platform.python_version(),sympy=s.__version__),
        input_sha256={'code/ordering_metric.py':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    target.write_text(json.dumps(payload,indent=2)+'\n')
    print('checks',len(tests),'phase norm',norm,c1,c2,'mass01',mat[0,1],flush=True)
if __name__=='__main__':main()
