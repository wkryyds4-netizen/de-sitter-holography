"""Merged transverse series and uniform-through-zero exact-arithmetic bounds.

The local distribution theorem uses Lebesgue dr, not a completed physical
mass integral. All evaluated tail bounds exclude numerical rounding.
"""
from fractions import Fraction
import sympy as sp
import mpmath as mp
from ds3_full_tail_bound import omega,numeric
from ds3_shadow_duality import log_reflection
from ds3_crossing_dictionary import period
from global_racah import shadow_coefficient
from ds3_transverse_phase import leading


def constants(M,u,v,p,q,R):
    M,u,v,p,q,R=map(mp.mpf,[M,u,v,p,q,R])
    alpha=[(v+q)/2,(p-u)/2,(-u-p)/2,(v-q)/2]
    if R<=0 or min(p,q)<=0 or max(abs(x)+R/2 for x in alpha)>mp.mpf('.5')+32*mp.eps:
        raise ValueError('Requires a closed transverse interval inside the numerator window')
    B=max(abs(u),abs(v))+R;m=M-B
    if m<=0:raise ValueError('Uniform lower mass must be positive')
    C0=mp.cosh(mp.pi/2)*mp.sqrt(1+R/(2*m))*mp.exp(mp.pi*R/4)/mp.sqrt(-mp.expm1(-2*mp.pi*m))
    P1=lambda t:mp.euler+1+mp.asinh(t)
    Ph=mp.euler+2*mp.log(2)+2+mp.asinh(1)
    Q=P1(R/2)+2*P1(M+B)+2*Ph
    RF=(mp.zeta(2)-1)/(1+m*m)
    RY=2*C0*((mp.zeta(2)-1)*Q+14)/(1+m*m)
    Smax=2*C0*Q;Ymax=Smax+RY
    Dmax=C0/(4*mp.pi)*mp.sqrt(mp.sinh(mp.pi*R/2)/(mp.pi*R/2))
    return dict(M=M,u=u,v=v,p=p,q=q,R=R,m=m,C0=C0,Q=Q,RF=RF,RY=RY,Smax=Smax,Ymax=Ymax,Dmax=Dmax)


def regularized_logs(M,u,v,p,q,r):
    M,u,v,p,q,r=map(mp.mpf,[M,u,v,p,q,r]);b=M+u;d0=M+v;z=-1j*r/2
    A=[(1+1j*x)/2 for x in [v+q,p-u,-u-p,v-q]]
    la=mp.loggamma(1+z)+mp.loggamma(1-1j*b)-mp.loggamma(1-1j*b-z)
    la+=sum(mp.loggamma(x-z)-mp.loggamma(x) for x in A[2:])
    lb=mp.loggamma(1-z)+mp.loggamma(1+1j*(d0+r))-mp.loggamma(1+1j*d0-z)
    lb+=sum(mp.loggamma(x)-mp.loggamma(x-z) for x in A[:2])
    return la,lb,A,z,b,d0


def merged_y(M,u,v,p,q,r,R,N=128):
    co=constants(M,u,v,p,q,R);r=mp.mpf(r)
    if abs(r)>co['R']:raise ValueError('Outside interval')
    la,lb,A,z,b,d0=regularized_logs(M,u,v,p,q,r)
    values=[];term=mp.mpc(1)
    if r==0:
        bracket=mp.digamma(1-1j*b)+mp.digamma(1+1j*d0)+2*mp.digamma(1)-sum(mp.digamma(a) for a in A)
        S=bracket
        for n in range(N):
            values.append(term*bracket)
            term*=mp.fprod(n+a for a in A)/((n+1+1j*d0)*(n+1-1j*b)*(n+1)**2)
            bracket+=1/(n+1-1j*b)+1/(n+1+1j*d0)+2/mp.mpf(n+1)-sum(1/(n+a) for a in A)
    else:
        # The two gamma poles cancel before summation. log-ratio branches
        # are continuous from r=0; expm1 preserves their small difference.
        lr=la-lb;B=mp.exp(lb);S=B*mp.expm1(lr)/z
        ba=[1+1j*(d0+r),1-1j*b-z,1-z]
        bb=[1+1j*d0-z,1-1j*b,1+z]
        for n in range(N):
            values.append(B*term*mp.expm1(lr)/z)
            lr+=sum(mp.log1p(-z/(n+a)) for a in A)
            lr+=sum(mp.log1p((y-x)/(n+x)) for x,y in zip(ba,bb))
            term*=mp.fprod(n+a for a in A)/(mp.fprod(n+x for x in bb)*(n+1))
    P=mp.fprod(mp.mpf(j*j)/(j*j+co['m']**2) for j in range(1,N+1))
    tail=2*co['C0']*P/N*(co['Q']+7*(2+mp.log(N)))
    return mp.fsum(values),tail,dict(S=S,N=N,constants=co,product=P,
                                    absolute_sum=mp.fsum(abs(x) for x in values),logs=[la,lb])


def evaluate(M,u,v,p,q,r,R='.25',N=128,refinement=2):
    vals=list(map(mp.mpf,[M,u,v,p,q,r]));M,u,v,p,q,r=vals;b=M+u;d=M+v+r
    hs=[(1+sp.I*sp.Rational(str(x)))/2 for x in [M+u+v,q,M,p]]
    h=(1+sp.I*sp.Rational(str(b)))/2;hp=(1+sp.I*sp.Rational(str(d)))/2
    X,ex,xd=omega(h,hp,hs,refinement)
    Y,ey,yd=merged_y(*vals,R,N)
    h,hp=map(numeric,[h,hp]);h1,h2,h3,h4=hn=list(map(numeric,hs))
    pref=mp.gamma(h+h1-h2)*mp.gamma(h-h1+h2)*mp.gamma(1-h+h4-h3)*mp.gamma(1-h+h3-h4)/(4*mp.gamma(2*h)*mp.gamma(2-2*h))
    cp=shadow_coefficient(1-hp,h1,h4)*pref*mp.sin(mp.pi*(hp-h1-h4))*mp.sin(mp.pi*(hp-h2-h3))
    cm=shadow_coefficient(hp,h3,h2)*pref*mp.sin(mp.pi*(1-hp-h1-h4))*mp.sin(mp.pi*(1-hp-h2-h3))
    cm=mp.exp(mp.log(cm)+log_reflection(h,1-hp,hn)+log_reflection(1-h,1-hp,hn))
    ns=period(M+u+v,q,b)*period(M,p,-b);nt=period(M,q,d)*period(M+u+v,p,-d)
    factor=b*d/mp.pi**4*(ns/abs(ns))/(nt/abs(nt));cp*=factor;cm*=factor
    value=cp*X*Y+cm*mp.conj(X*Y)
    tail=(abs(cp)+abs(cm))*(abs(X)*ey+abs(Y)*ex+ex*ey)
    ca,cb=xd['coefficients'];D0=cp*ca;E0=cp*cb;co=yd['constants']
    alpha=[(v+q)/2,(p-u)/2,(-u-p)/2,(v-q)/2];B=b-abs(r)/2;astar=max(map(abs,alpha))+abs(r)/2
    Bb=mp.sqrt(mp.sinh(mp.pi*B)/(mp.pi*B))*mp.exp(2*mp.pi*astar+mp.pi*astar*astar/B)
    poly=list(map(Fraction,xd['groups'][1]['certificate']['coefficients']));K=max(abs(x) for x in poly[1:]);Nstar=2+(K.numerator+K.denominator-1)//K.denominator
    BG=(3*Nstar+3)*Bb
    bound=2*co['Dmax']*(co['RY']+co['RF']*co['Ymax'])+2*abs(E0)*BG*co['Ymax']
    finite=2*mp.re(D0*yd['S']);asym=leading(*vals)
    return dict(M=M,u=u,v=v,p=p,q=q,r=r,R=mp.mpf(R),value=value,tail=tail,X=X,Y=Y,
                Xtail=ex,Ytail=ey,D0=D0,E0=E0,coefficients=[cp,cm],finite=finite,finite_error=abs(value-finite),
                uniform_finite_bound=bound,passed=abs(value-finite)<=bound+tail,
                sine=asym,sine_error=abs(value-asym),scaled_error=abs(value-asym)*M/mp.log(M),
                BG=BG,Bb=Bb,Nstar=Nstar,Xdetails=xd,Ydetails=yd)
