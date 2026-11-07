"""Exact finite Racah moments and a symbolic uniform sixth-moment certificate.

The symbolic expansion uses e=N^(-1/2), a=A/e, b=B/e,
m=1/(2e^2)+d. Coefficients lie in QQ[A,B,d], not sampled parameters.
"""
from pathlib import Path
from fractions import Fraction as F
from math import factorial, comb
import hashlib,json,platform,time
import sympy as sp
from sympy.polys.rings import ring
from sympy.polys.domains import QQ
from cubic_linear_dynamics import crossing

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-filter-limit'


def finite_deficit_moments(N,m,a,b,order=6):
    assert order<=min(m,N-m)
    x=lambda j:j*(N+1-j)
    polys=[[1]];fall=[F(1)];raw=[F(1)]
    for n in range(1,order+1):
        old=polys[-1];new=[0]*(len(old)+1)
        for j,c in enumerate(old):new[j]-=x(n-1)*c;new[j+1]+=c
        polys.append(new)
        coeff=[F(1)]
        for j in range(n):
            coeff.append(coeff[-1]*F((-n+j)*(n-N-1+j),(-m+j)*(m-N+j)*(-N-1+j)*(j+1)))
        value=(crossing(N,m,a,n)*crossing(N,m,b,n)-sum(coeff[j]*fall[j] for j in range(n)))/coeff[n]
        fall.append(value)
        raw.append(value-sum(F(polys[n][j])*raw[j] for j in range(n)))
    top=x(a+b)
    return [sum(F(comb(k,j)*(-1)**j*top**(k-j))*raw[j] for j in range(k+1))/N**k for k in range(order+1)]


def symbolic_certificate():
    R,A,B,d=ring('A,B,d',QQ);zero,one=R.zero,R.one;T=12
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
            for v in [2*(d-k),-2*(d+k),1-k]:den=mul(den,quad(v))
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
        z=sp.Symbol('z');touch=sp.Integer(1)
        for _ in range(n):touch=sp.expand(z*(touch+sp.diff(touch,z)))
        expected=R.from_expr(touch.subs(z,2*sp.Symbol('A')*sp.Symbol('B')))
        assert delta[n]==expected,(n,delta[n],expected)
        rows.append(dict(order=n,cancelled_negative_powers=n,limit=str(delta[n]),
                         coefficient_parameter_domain='QQ[A,B,d]',charge_shift_d_cancels=True))
        print('symbolic uniform moment',n,str(delta[n]),flush=True)
    return rows


def main():
    OUT.mkdir(exist_ok=True);target=OUT/'sixth-moment-certificate.json'
    if target.exists():raise FileExistsError(target)
    start=time.time();rows=symbolic_certificate()
    names=['code/memory_sixth_moment.py','code/cubic_linear_dynamics.py']
    result=dict(status='passed',symbolic_rows=rows,formal_series_order=12,
                total_cancelled_negative_powers=sum(r['cancelled_negative_powers'] for r in rows),
                scope='Exact polynomial-ring identities for all A,B,d. Analytic rational denominators stay nonzero near e=0 on compact A,B,d sets; hence scaled moments through sixth are uniformly bounded. Does not assert a full Poisson law.',
                elapsed_seconds=time.time()-start,environment=dict(python=platform.python_version(),sympy=sp.__version__),
                input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names})
    target.write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
