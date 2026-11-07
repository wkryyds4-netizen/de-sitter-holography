"""Independent exact integer-complex CAR check and saved symbolic certificates."""
import json,hashlib
from itertools import combinations
import sympy as s
from cubic_singlet import ROOT,OUT,q,z

def main():
    target=OUT/'identity-check.json'
    if target.exists():raise FileExistsError(target)
    checks=[]
    for n,p in [(3,1),(4,2)]:
        d=1<<n;cs=[]
        for i in range(n):
            c=s.zeros(d)
            for a in range(d):
                if a>>i&1:c[a^(1<<i),a]=(-1)**((a&((1<<i)-1)).bit_count())
            cs.append(c)
        basis=[]
        for inds in combinations(range(n),p):
            b=s.eye(d)
            for i in inds:b=b*cs[i].T
            basis.append(b)
        H=s.zeros(d)
        for i,a in enumerate(basis):
            for j,b in enumerate(basis):H+=((i+j)%3-1+s.I*(i-j))*a*b.T
        L=lambda a:H*a-a*H
        fs=[[c,L(c),L(L(c)),L(L(L(c)))] for c in cs]
        B=lambda a,b:sum((v[a].conjugate().T*v[b] for v in fs),s.zeros(d))
        Q=s.diag(*[a.bit_count() for a in range(d)])
        D=lambda a:Q*a-sum((c.T*a*c for c in cs),s.zeros(d))
        A=B(1,1);M=B(0,3);T=(M+M.conjugate().T)/2
        tests={'M1':B(0,1)+p*H,'M2':B(0,2)-A,'M3':M-2*B(1,2)+B(2,1),
         'anti_Hermitian':M-M.conjugate().T-3*L(A),'rank_A':A-2*p*H**2+D(H**2),
         'rank_T':T-3*p*H**3+D(H**3)+s.Rational(3,2)*(H*A+A*H)}
        for name,val in tests.items():
            assert val.applyfunc(s.expand)==s.zeros(d),(n,p,name)
            checks.append(dict(N=n,p=p,identity=name,passed=True))
        print('Exact integer-complex CAR',n,p,'passed',flush=True)
    data=json.loads((OUT/'study.json').read_text());parse=lambda v:s.sympify(v.replace('lambda','lam'),locals={'q':q,'lam':z})
    G=s.Matrix([[parse(x) for x in row] for row in data['gram']]);K=s.Matrix([[parse(x) for x in row] for row in data['kinetic']])
    v=s.Matrix([1]+[-parse(c) for c in data['projection']]);V=parse(data['residual_variance'])
    B=1-q**3+4*z*q**2*(2*q**2-q+2)
    algebra={
      'projection_S':(v.T*G)[1], 'projection_Z':(v.T*G)[2],
      'positive_variance':V-4*z*(1+q)*B/((1+q*q)*(1+q+q*q)),
      'kinetic_ratio':(v.T*K*v)[0]-2*(1-q**3)*V,
      'small_lambda_variance':s.limit(V.subs(q,s.exp(-4*z))/z**2,z,0)-32,
      'zero_kinetic_Z':sum(K[i,2]**2 for i in range(3))}
    for name,val in algebra.items():assert s.simplify(val)==0,name
    spec=json.loads((OUT/'spectrum-check.json').read_text())
    assert all(r['matches'] and r['q3_difference']=='0' for r in spec['rows'])
    assert all(r['R_S']==r['R_Z']=='0' for r in spec['cross_moments'])
    target.write_text(json.dumps(dict(status='passed',exact_CAR=checks,symbolic_certificates=list(algebra),
     two_arc_matches=len(spec['rows']),vanishing_cross_moments=2*len(spec['cross_moments']),
     environment=dict(sympy=s.__version__),
     input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in ['code/check_cubic_identities.py','results/susskind-cubic-singlet/study.json','results/susskind-cubic-singlet/spectrum-check.json']}),indent=2)+'\n')

if __name__=='__main__':main()
