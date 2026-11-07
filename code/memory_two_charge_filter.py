"""Candidate exact two-charge filters derived from charged low/high eigenvectors."""
from fractions import Fraction as F
from math import comb

def choose(n,k):return comb(n,k) if 0<=k<=n else 0

def filters(N,p,m,a,b,r):
    assert m>p and 0<=a<=b<=min(p,m,N-m)
    d=m-p;n=N-m+1
    C=choose(m,p)*choose(N-m+p,p)
    Cm=choose(m-1,p)*choose(N-m+p+1,p)
    jc=choose(m-1,p)*choose(N-m+p,p)
    def kap(j,charge):return choose(charge-j,p)*choose(N-charge+p-j,p)
    def channel(j):
        low=choose(m-j,p)*choose(N-m+p-j+1,p)
        high=choose(m-j-1,p)*choose(N-m+p-j,p)
        e=F(low-high,N+1-2*j)
        alpha=high+e*(n-j)
        beta=-e*F((n-j)*(m-j),d)
        gamma=-d*e
        delta=high+e*(m-j)
        u=m*gamma+d*(delta+jc)
        z=F((n-j)*(m-j),d)
        v=high-j*e
        return dict(alpha=alpha,beta=beta,gamma=gamma,delta=delta,u=u,z=z,v=v,e=e)
    aa,bb,rr=channel(a),channel(b),channel(r)
    D4=(F(0),F(3*Cm+kap(a,m-1)+kap(b,m-1)+kap(r,m-1)))
    D13=(d*(4*Cm+kap(a,m-1)+kap(b,m-1)+aa['gamma']+bb['gamma']),
         2*(jc+rr['v'])+aa['delta']+bb['delta'])
    D22=(2*m*Cm+aa['u']+bb['u']+d*(aa['alpha']+bb['alpha'])+2*d*rr['alpha'],
         2*C+aa['beta']+bb['beta']-rr['e']*(aa['z']+bb['z']))
    D121=(aa['u']+bb['u']+m*Cm,F(kap(r,m)))
    D40=(F(3*C+kap(a,m)+kap(b,m)+kap(r,m)),F(0))
    words=dict(D4=D4,D13=D13,D22=D22,D121=D121,D40=D40)
    g=tuple(D4[i]-2*D13[i]+F(5,4)*D22[i]+F(7,2)*D121[i]+(4*p-3*m)*D40[i] for i in range(2))
    return dict(g=g,words=words,channels=[aa,bb,rr])
