"""Termwise merged LPRS Gamma-pair series including its exact collision.

All scalar coefficients use arbitrary-precision division. The first
collision prototype is retained separately as historical evidence.
"""
import mpmath as mp
from ds3_racah_collision import collision_data
from global_racah import shadow_coefficient
from ds3_crossing_dictionary import period


def asymptotic(top,bottom,K):
    logs=[mp.mpc(0)];cs=[mp.mpc(1)]
    for k in range(1,K+1):
        diff=sum(mp.bernpoly(k+1,a) for a in top)-sum(mp.bernpoly(k+1,b) for b in bottom)
        logs.append((-1)**(k+1)*diff/(k*(k+1)))
        cs.append(sum(j*logs[j]*cs[k-j] for j in range(1,k+1))/k)
    return cs


def omega_merged(h,hp,hs,N=128,K=24):
    hp0,top,bottom=collision_data(h,hs);z=hp0-hp
    hit=abs(z)<64*mp.eps*max(1,abs(h),abs(hp),*[abs(x) for x in hs])
    values=[]
    if hit:
        term=mp.mpc(1)
        for n in range(N):
            bracket=mp.digamma(2*h+n)+mp.digamma(2*hp0+n)+2*mp.digamma(1+n)-sum(mp.digamma(a+n) for a in top)
            values.append(term*bracket)
            term*=mp.fprod(a+n for a in top)/(mp.fprod(b+n for b in bottom)*(n+1))
        cs=asymptotic(top,bottom+[mp.mpf(1)],K)
        ds=[mp.mpc(0)]
        for k in range(1,K+1):
            diff=mp.bernpoly(k,2*h)+mp.bernpoly(k,2*hp0)+2*mp.bernpoly(k,1)-sum(mp.bernpoly(k,a) for a in top)
            ds.append((-1)**(k+1)*diff/k)
        coeff=[sum(cs[j]*ds[k-j] for j in range(k)) for k in range(1,K+1)]
        C=mp.fprod(mp.gamma(b) for b in bottom)/mp.fprod(mp.gamma(a) for a in top)
    else:
        ta=[a-z for a in top];ba=[2*hp0-2*z,2*h-z,1-z,mp.mpf(1)]
        tb=top;bb=[2*hp0-z,2*h,1+z,mp.mpf(1)]
        C=mp.gamma(1+z)*mp.gamma(1-z)*mp.gamma(2*hp0-2*z)*mp.gamma(2*h)/(mp.gamma(top[0]-z)*mp.gamma(top[1]-z)*mp.gamma(top[2])*mp.gamma(top[3]))
        B=mp.fprod(mp.gamma(a) for a in tb)/mp.fprod(mp.gamma(b) for b in bb)
        lr=sum(mp.loggamma(a) for a in ta)-sum(mp.loggamma(b) for b in ba)-sum(mp.loggamma(a) for a in tb)+sum(mp.loggamma(b) for b in bb)
        # Reduce only a constant 2*pi*i branch multiple, preserving exp(lr).
        lr-=2j*mp.pi*mp.nint(mp.im(lr)/(2*mp.pi))
        for n in range(N):
            values.append(C*B*mp.expm1(lr)/z)
            lr+=sum(mp.log1p(-z/(n+a)) for a in top)
            lr+=mp.log1p(z/(n+2*hp0-2*z))-mp.log1p(-z/(n+2*h))+mp.log1p(2*z/(n+1-z))
            B*=mp.fprod(n+a for a in tb)/mp.fprod(n+b for b in bb)
        ca=asymptotic(ta,ba,K);cb=asymptotic(tb,bb,K)
        coeff=[(ca[k]-cb[k])/z for k in range(1,K+1)]
    tails=[C*v*mp.zeta(k+2,N) for k,v in enumerate(coeff,1)]
    value=mp.fsum(values)+mp.fsum(tails)
    return value,dict(collision=hit,z=z,N=N,K=K,tail=mp.fsum(tails),last_tail_term=tails[-1],
                     summation_condition=(mp.fsum(abs(v) for v in values)+mp.fsum(abs(v) for v in tails))/abs(value),
                     leading_tail_coefficient=coeff[0])


def physical_merged(b,d,external,N=128,K=24):
    b,d=mp.mpf(b),mp.mpf(d);external=list(map(mp.mpf,external))
    hs=[(1+1j*x)/2 for x in external];h,hp=(1+1j*b)/2,(1+1j*d)/2;h1,h2,h3,h4=hs
    pref=mp.gamma(h+h1-h2)*mp.gamma(h-h1+h2)*mp.gamma(1-h+h4-h3)*mp.gamma(1-h+h3-h4)/(4*mp.gamma(2*h)*mp.gamma(2-2*h))
    terms=[];records=[]
    for q,co in [(hp,shadow_coefficient(1-hp,h1,h4)),(1-hp,shadow_coefficient(hp,h3,h2))]:
        vs=[]
        for H in [h,1-h]:
            v,r=omega_merged(H,q,hs,N,K);vs.append(v);records.append(r)
        terms.append(co*pref*mp.sin(mp.pi*(q-h1-h4))*mp.sin(mp.pi*(q-h2-h3))*vs[0]*vs[1])
    a,q,c,p=external;ns=period(a,q,b)*period(c,p,-b);nt=period(c,q,d)*period(a,p,-d)
    value=b*d/mp.pi**4*(ns/abs(ns))/(nt/abs(nt))*mp.fsum(terms)
    return value,dict(omega=records,shadow_condition=mp.fsum(abs(v) for v in terms)/abs(mp.fsum(terms)))
