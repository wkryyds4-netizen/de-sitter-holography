"""Central-window Racah certificate: m=N/2+z sqrt(N), |z| bounded.
Adapted from the bounded-charge-shift recurrence without changing old inputs.
"""
from pathlib import Path
from fractions import Fraction as F
from math import factorial, comb
import hashlib,json,platform,time
import sympy as sp
from sympy.polys.rings import ring
from sympy.polys.domains import QQ
from memory_sixth_moment import finite_deficit_moments
from rank_product_moments import product_distribution
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-central-window'

def symbolic_certificate():
    R,A,B,z=ring('A,B,z',QQ);zero,one=R.zero,R.one;T=12
    def const(c):return [R(c)]+[zero]*T
    def trim(x):return (x+[zero]*(T+1))[:T+1]
    def add(x,y):return [a+b for a,b in zip(trim(x),trim(y))]
    def neg(x):return [-v for v in x]
    def sub(x,y):return add(x,neg(y))
    def mul(x,y):
        out=[zero]*(T+1)
        for i,a in enumerate(x):
            if a:
                for j,b in enumerate(y[:T+1-i]):
                    if b:out[i+j]+=a*b
        return out
    def power(x,n):
        out=const(1)
        for _ in range(n):out=mul(out,x)
        return out
    def div(x,y,limit=T):
        c=y[0].get((0,0,0),QQ.zero)
        assert c and y[0]==R(c)
        out=[zero]*(T+1)
        for i in range(limit+1):
            out[i]=(x[i]-sum((y[j]*out[i-j] for j in range(1,i+1)),zero)).quo_ground(c)
        return out
    def quad(c):return [one,zero,R(c)]+[zero]*(T-2)
    def coefficient(n,j):
        num=const(1);den=const(1)
        for k in range(j):
            num=mul(num,quad(-(n-1+k)))
            den=mul(den,trim([one,2*z,R(-2*k)]))
            den=mul(den,trim([one,-2*z,R(-2*k)]))
            den=mul(den,quad(1-k))
        lead=F((-1)**j*factorial(n),factorial(n-j)*factorial(j))*4**j
        ratio=div(num,den)
        return [zero]*j+[v*QQ(lead.numerator,lead.denominator) for v in ratio[:T+1-j]]
    roots=[trim([zero,R(l),zero,R(l*(1-l))]) for l in range(6)]
    # Polynomials in Y, with power-series coefficients.
    polys=[[const(1)]]
    for n in range(1,7):
        old=polys[-1];new=[const(0) for _ in range(n+1)]
        for j,c in enumerate(old):new[j]=sub(new[j],mul(roots[n-1],c));new[j+1]=add(new[j+1],c)
        polys.append(new)
    ya=trim([A,-A*A,A]);yb=trim([B,-B*B,B]);yt=trim([A+B,-(A+B)**2,A+B])
    def at_y(n,y):
        out=const(1)
        for j in range(n):out=mul(out,sub(y,roots[j]))
        return out
    fall=[const(1)];raw=[const(1)];rows=[]
    for n in range(1,7):
        cs=[coefficient(n,j) for j in range(n+1)]
        ca=const(0);cb=const(0)
        for j in range(n+1):ca=add(ca,mul(cs[j],at_y(j,ya)));cb=add(cb,mul(cs[j],at_y(j,yb)))
        numerator=mul(ca,cb)
        for j in range(n):numerator=sub(numerator,mul(cs[j],fall[j]))
        assert all(v==zero for v in numerator[:n]),('triangular valuation',n)
        value=div(trim(numerator[n:]),trim(cs[n][n:]),T-n)
        fall.append(value)
        moment=value
        for j in range(n):moment=sub(moment,mul(polys[n][j],raw[j]))
        raw.append(moment[:T-n+1]+[zero]*n)
        delta=const(0)
        for j in range(n+1):delta=add(delta,[v*((-1)**j*comb(n,j)) for v in mul(power(yt,n-j),raw[j])])
        assert all(v==zero for v in delta[:n]),('negative power',n,[str(v) for v in delta[:n]])
        tvar=sp.Symbol('t');touch=sp.Integer(1)
        for _ in range(n):touch=sp.expand(tvar*(touch+sp.diff(touch,tvar)))
        expected=R.from_expr(touch.subs(tvar,2*sp.Symbol('A')*sp.Symbol('B')))
        assert delta[n]==expected,(n,delta[n],expected)
        rows.append(dict(order=n,cancelled_negative_powers=n,limit=str(delta[n]),
                         coefficient_parameter_domain='QQ[A,B,z]',central_charge_z_cancels=True))
        print('central uniform moment',n,str(delta[n]),flush=True)
    return rows


def main():
    OUT.mkdir(exist_ok=True)
    target=OUT/'central-moments.json'
    if target.exists():raise FileExistsError(target)
    start=time.time(); rows=symbolic_certificate()
    small=[]
    for N,m,a,b in [(16,6,2,3),(16,10,2,3),(20,7,3,3),(20,13,3,3)]:
        pi=product_distribution(N,m,a,b)
        moment=finite_deficit_moments(N,m,a,b)
        x=lambda j:j*(N+1-j)
        direct=[sum(prob*F(x(a+b)-x(r),N)**k for r,prob in enumerate(pi)) for k in range(7)]
        assert moment==direct
        sixthS=sum(prob*(a+b-r)**6 for r,prob in enumerate(pi))
        bound=moment[6]*F(N,N+1-2*(a+b))**6
        assert sixthS<=bound
        small.append(dict(N=N,m=m,a=a,b=b,moments=list(map(str,moment)),sixth_S=str(sixthS),bound=str(bound)))
    finite=[]
    for p in [16,32,64,128]:
        N=4*p*p
        for z in [-2,-1,0,1,2]:
            m=N//2+2*p*z
            for u,v in [(0,0),(1,2)]:
                values=finite_deficit_moments(N,m,p-u,p-v)
                opposite=finite_deficit_moments(N,N-m,p-u,p-v)
                assert values==opposite
                t=F(2*(p-u)*(p-v),N)
                target6=t**6+15*t**5+65*t**4+90*t**3+31*t*t+t
                bound=values[6]*F(N,N+1-4*p)**6
                finite.append(dict(N=N,p=p,z=z,input_deficits=[u,v],sixth=str(values[6]),
                    sixth_decimal=float(values[6]),finite_scaled_touchard=float(target6),
                    deviation=float(values[6]-target6),S_sixth_upper_bound=str(bound)))
    names=['code/memory_central_window.py','code/memory_sixth_moment.py',
           'code/cubic_linear_dynamics.py','code/rank_product_moments.py']
    data=dict(status='passed',symbolic_rows=rows,series_order=12,
        exact_full_distribution_checks=small,finite_cases=finite,
        cancelled_negative_powers=sum(r['cancelled_negative_powers'] for r in rows),
        scope='Uniform sixth moment for independent isotropic rank products on bounded z=(m-N/2)/sqrt(N) windows and bounded a/sqrt(N),b/sqrt(N). Not the complete memory Gram or a growing-z window.',
        environment=dict(python=platform.python_version(),sympy=sp.__version__),elapsed_seconds=time.time()-start,
        input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names})
    target.write_text(json.dumps(data,indent=2)+'\n')
    print('passed central certificate and',len(finite),'finite cases',flush=True)

if __name__=='__main__':main()
