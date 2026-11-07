"""Exact signed rank-product moments at half filling and their Gram consequences."""
from pathlib import Path
from fractions import Fraction as F
from math import comb
import json,hashlib,platform
import sympy as s
from cubic_linear_dynamics import crossing

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-parity-gram'

def rates(N,j):
    m=N//2
    B=F((j-m)*(j-N-1)**2,2*(2*j-N-1))
    D=F(j*j*(j-m-1),2*(2*j-N-1))
    return B,D

def product_table(N):
    m=N//2;dim=comb(N,m)
    ds=[comb(N,j)**2-(comb(N,j-1)**2 if j else 0) for j in range(m+1)]
    chi=[[crossing(N,m,a,b) for b in range(m+1)] for a in range(m+1)]
    pi={}
    for a in range(m+1):
        for b in range(m+1):
            values=[F(ds[r],dim*dim)*sum(ds[j]*chi[a][j]*chi[b][j]*chi[r][j] for j in range(m+1)) for r in range(m+1)]
            assert min(values)>=0 and sum(values)==1
            pi[a,b]=values
    return ds,chi,pi

def signed_dx(N,a,b):
    m=N//2;x=lambda j:j*(N+1-j)
    f=crossing(N,m,a,b)
    fm=crossing(N,m,a-1,b) if a else f
    return [f,-2*a*b*f,(4*a*a*b*b-2*a*x(b))*f+a*a*(m+1-a)*(fm-f)]

def main():
    OUT.mkdir(exist_ok=True);target=OUT/'signed-moments.json'
    if target.exists():raise FileExistsError(target)
    rows=[];moments=0;recurrences=0
    for N in [2,4,6,8,10,12]:
        m=N//2;dim=comb(N,m);ds,chi,pi=product_table(N);x=lambda j:j*(N+1-j)
        for b in range(m+1):
            assert sum((-1)**a*ds[a]*chi[a][b] for a in range(m+1))==(-1)**(m+b)*dim
            for a in range(m+1):
                B,D=rates(N,a)
                L=B*(chi[a][b]-(chi[a+1][b] if a<m else 0))+D*(chi[a][b]-(chi[a-1][b] if a else 0))
                assert L==x(b)*chi[a][b]
                assert 2*(B+D)==F(N*(N+1),2)-x(a)
                recurrences+=1
                actual=[sum((-1 if (a+b-r)%2 else 1)*pi[a,b][r]*(x(a+b)-x(r))**k for r in range(m+1)) for k in range(3)]
                assert actual==signed_dx(N,a,b),(N,a,b,actual,signed_dx(N,a,b))
                moments+=3
        rows.append(dict(N=N,m=m,input_rank_pairs=(m+1)**2,signed_moment_checks=3*(m+1)**2))
    # Symbolic operator manipulation behind the stable signed second moment.
    N,a,b=s.symbols('N a b');m=N/2
    x=lambda j:j*(N+1-j)
    K=N*(N+1)/2;t=N-2*a
    D=a*a*(a-m-1)/(2*(2*a-N-1))
    identity=-t*(K-x(a))/2+t*x(b)+N*(K/2-x(a)-x(b))+2*D*(t+1)
    assert s.factor(identity+2*a*x(b))==0
    lam,q=s.symbols('lam q',positive=True);theta=2*lam;z=1+q
    G=s.Matrix([[z,0,theta*(1-q)],[0,2*lam*z,0],
        [theta*(1-q),0,theta**2*z+theta*(1-q)]])
    residual=s.factor(G[2,2]-G[0,2]**2/G[0,0])
    assert s.factor(residual-(theta*(1-q)+4*theta**2*q/z))==0
    small=s.limit((residual/lam**2).subs(q,s.exp(-4*lam)),lam,0)
    assert small==16
    names=['code/memory_parity_moments.py','code/cubic_linear_dynamics.py']
    result=dict(status='passed',rows=rows,signed_moment_checks=moments,
        recurrence_checks=recurrences,
        exact_signed_delta_moments=['chi(a,b)','-2*a*b*chi(a,b)',
            '(4*a^2*b^2-2*a*x_b)*chi(a,b)+a^2*(m+1-a)*(chi(a-1,b)-chi(a,b))'],
        gram_over_C_squared=[[str(x) for x in row] for row in G.tolist()],
        residual_over_C_squared=str(residual),small_lambda_residual_over_Js_four=str(small),
        scope='Exact half filling; independent ordered-product law. Gram is the same-coupling symmetric second-chaos limit, not a time-dependent spectrum.',
        environment=dict(python=platform.python_version(),sympy=s.__version__),
        input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names})
    target.write_text(json.dumps(result,indent=2)+'\n')
    print('passed',moments,'signed moments,',recurrences,'recurrences; small lambda',small)

if __name__=='__main__':main()
