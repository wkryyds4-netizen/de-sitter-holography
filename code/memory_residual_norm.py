"""Exact half-filled second-chaos A norm, retaining adjacent-charge coherence."""
from fractions import Fraction as F
from cubic_linear_dynamics import crossing, weights
from memory_coherent_cross import exchanged


def pair_norm(N,p,a,b):
    """Contribution per ordered rank pair, divided by C_a C_b."""
    assert N%2==0
    m=N//2;n=m+1
    assert p<m and a<=p and b<=p
    xa=a*(N+1-a);xb=b*(N+1-b)
    g=F((m-p)**2,m);h=2*p-m;eta=F(p*p,m)
    P=g*g*F(xa*xb,(m*n-xa)*(m*n-xb))
    cf=F((m-1)*n*(N+1),N)
    chi=crossing(N,m,a,b);delta=exchanged(N,a,b)-chi
    value=eta*eta*(1+chi)+P*(F(m*n,1)/cf-1-chi)+(g*g/n-P/n+2*h*g)*delta
    # A norm is positive even though the exact summands need not be.
    assert value>=0,(N,p,a,b,value)
    return value


def pair_norm_raw(N,p,a,b):
    m=N//2;n=m+1
    x=lambda r:r*(N+1-r)
    lower_cf=F((m-1)*n*(N+1),N)
    ratio=lambda r:F((m-p)**2*n,m*(m*n-x(r)))
    chiminus=crossing(N,m-1,a,b)
    rem=ratio(a)*ratio(b)*F(m,n)*((m*n-x(a)-x(b))*(1+chiminus)+F(x(a)*x(b),1)/lower_cf)
    g=F((m-p)**2,m);h=2*p-m
    return rem+2*h*g*(1+exchanged(N,a,b))+h*h*(1+crossing(N,m,a,b))


def exact_norm(N,p):
    omega=weights(N,N//2,p)
    return sum(omega[a]*omega[b]*pair_norm(N,p,a,b) for a in range(p+1) for b in range(p+1))
