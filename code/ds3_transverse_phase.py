"""Full-kernel transverse leading term in the established G=0 CG convention.

Exact-arithmetic bounds only. Punctured-window asymptotics and fixed-lambda
scaling do not establish a distributional or physical wave-packet limit.
"""
from fractions import Fraction
import mpmath as mp
from ds3_shadow_duality import physical
from ds3_controlled_omega import controlled


def scattering_phase(u, v, p, q, r):
    u, v, p, q, r = map(mp.mpf, [u, v, p, q, r])
    A = [(1+1j*x)/2 for x in [v+q, p-u, -u-p, v-q]]
    z = -1j*r/2
    # An unwrapped real phase, continuous from delta(0)=0.
    delta = mp.im(2*mp.loggamma(1+z)
                  + sum(mp.loggamma(a-z)-mp.loggamma(a) for a in A))
    return delta


def leading(M, u, v, p, q, r):
    M, u, v, p, q, r = map(mp.mpf, [M, u, v, p, q, r])
    if r == 0:
        A = [(1+1j*x)/2 for x in [v+q, p-u, -u-p, v-q]]
        return (mp.log(M)+(2*mp.digamma(1)-mp.re(sum(mp.digamma(a) for a in A)))/2)/mp.pi
    return mp.sin(r*mp.log(M)-scattering_phase(u,v,p,q,r))/(mp.pi*r)


def evaluate(M, u, v, p, q, r, refinement=2):
    M,u,v,p,q,r=map(mp.mpf,[M,u,v,p,q,r])
    if r == 0:
        raise ValueError('Use the previously proved collision formula at r=0')
    b=M+u;d0=M+v;d=d0+r
    alpha=[(v+q)/2,(p-u)/2,(-u-p)/2,(v-q)/2]
    if max(abs(x) for x in alpha)>mp.mpf('.5') or max(abs(x+r/2) for x in alpha)>mp.mpf('.5'):
        raise ValueError('Outside the proved numerator window')
    value,tail,detail=physical(str(b),str(d),[str(M+u+v),str(q),str(M),str(p)],refinement)
    C=detail['coefficients'][0]
    xdata=detail['omega'][0]['details'];ydata=detail['omega'][1]['details']
    ca,cb=xdata['coefficients'];ya,yb=ydata['coefficients']
    D0=C*ca;E0=C*cb;S=ya+yb
    # Both F series use the small numerators alpha+r/2.
    ma=min(d,b-r/2);mb=min(d0+r/2,b)
    if min(ma,mb)<=0:raise ValueError('Mass denominator not positive')
    RF=(mp.zeta(2)-1)/(1+ma*ma)
    RY=(mp.zeta(2)-1)*(abs(ya)/(1+ma*ma)+abs(yb)/(1+mb*mb))
    # Whole companion G: a bound for every coefficient and a polynomial tail.
    B=b-abs(r)/2;Astar=max(abs(x) for x in alpha)+abs(r)/2
    assert B>0 and b+d0+r/2>=B
    Bb=mp.sqrt(mp.sinh(mp.pi*B)/(mp.pi*B))*mp.exp(2*mp.pi*Astar+mp.pi*Astar*Astar/B)
    polynomial=[Fraction(x) for x in xdata['groups'][1]['certificate']['coefficients']]
    assert polynomial[0]==1 and len(polynomial)==11
    K=max(abs(x) for x in polynomial[1:]);Nstar=2+(K.numerator+K.denominator-1)//K.denominator
    BG=(3*Nstar+3)*Bb
    finite=2*mp.re(D0*S)
    bound=2*abs(D0)*(RY+RF*(abs(S)+RY))+2*abs(E0)*BG*(abs(S)+RY)
    asym=leading(M,u,v,p,q,r)
    z=-1j*r/2;A=[mp.mpf('.5')+1j*x for x in alpha]
    Bsmall=mp.gamma(A[0])*mp.gamma(A[1])/(mp.gamma(A[0]-z)*mp.gamma(A[1]-z))
    amp=mp.exp(mp.pi*r/4)/(4*mp.pi*abs(mp.gamma(1-z)))
    amp*=mp.sqrt(mp.fprod(mp.cosh(mp.pi*(x+r/2))/mp.cosh(mp.pi*x) for x in alpha[2:]))
    phase=-mp.im(mp.loggamma(A[0])+mp.loggamma(A[1])-mp.loggamma(A[0]-z)-mp.loggamma(A[1]-z))+mp.im(mp.loggamma(1+z))
    Dinf=amp*mp.exp(1j*phase)
    Asmall=mp.gamma(A[2]-z)*mp.gamma(A[3]-z)/(mp.gamma(A[2])*mp.gamma(A[3]))
    T=Dinf*mp.gamma(z)*mp.exp(-1j*mp.pi*z/2)*Asmall
    J=Dinf*mp.gamma(-z)*mp.exp(-1j*mp.pi*z/2)*Bsmall
    expectedT=1j*mp.exp(1j*scattering_phase(u,v,p,q,r))/(2*mp.pi*r)
    return dict(M=M,u=u,v=v,p=p,q=q,r=r,value=value,tail=tail,
                finite_leading=finite,finite_error=abs(value-finite),finite_bound=bound,
                asymptotic_leading=asym,asymptotic_error=abs(value-asym),
                delta=scattering_phase(u,v,p,q,r),D0=D0,E0=E0,D_infinity=Dinf,
                D_scaled=D0*mp.exp(-z*mp.log(M)),T=T,J=J,
                scattering_identity_error=abs(T-expectedT),nonoscillatory_real=mp.re(J),
                RF=RF,RY=RY,S=S,BG=BG,Bb=Bb,Nstar=Nstar,
                finite_passed=abs(value-finite)<=bound+tail,
                coefficients=[ca,cb,ya,yb],details=detail)
