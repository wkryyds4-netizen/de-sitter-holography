"""General-charge Casimir-deficit cross with the coherent force."""
from fractions import Fraction as F
from cubic_linear_dynamics import crossing
from memory_central_norm import exchanged


def pair_cross(N,p,m,a,b):
    assert p<min(m,N-m+1) and max(a,b)<=p
    xa=a*(N+1-a);xb=b*(N+1-b)
    cf=F(m*(N-m)*(N+1),N)
    ct=F((m-1)*(N-m)*(N+1),N-1)
    chi=crossing(N,m,a,b);xi=exchanged(N,m,a,b);D=xi-chi
    g=F((m-p)**2,m);h=2*p-m;eta=F(p*p,m)
    neutral=-2*a*b*(1+chi)+F(xa*xb,cf)
    coherent=-2*a*b*(1+xi)+F(xa*xb,ct)+m*D
    raw=(g*coherent+h*neutral)/N
    stable=(eta*neutral+g*(F(xa*xb,m*(m-1)*(N+1))+(m-2*a*b)*D))/N
    assert raw==stable
    return neutral/N,stable
