"""General-charge coherent A variance and B cross; exact rational kernels."""
from fractions import Fraction as F
from cubic_linear_dynamics import crossing,weights


def exchanged(N,m,a,b):
    term=F(1);result=term
    for j in range(min(a,b)):
        term*=F((-a+j)*(a-N-1+j)*(-b+j)*(b-N-1+j),
                (1-m+j)*(m-N+j)*(-N-1+j)*(j+1))
        result+=term
    return result


def pair(N,p,m,a,b):
    assert p<min(m,N-m+1) and max(a,b)<=p
    n=N-m+1;xa=a*(N+1-a);xb=b*(N+1-b)
    g=F((m-p)**2,m);h=2*p-m;eta=F(p*p,m)
    P=g*g*F(xa*xb,(m*n-xa)*(m*n-xb))
    cf=F((m-1)*n*(N+1),N)
    chi=crossing(N,m,a,b);xi=exchanged(N,m,a,b);D=xi-chi
    lower=crossing(N,m-1,a,b)
    assert lower==chi+F(n-m,n)*D
    norm=eta*eta*(1+chi)+P*(F(m*n,1)/cf-1-chi)+((g*g-P)*F(n-m,n)+2*h*g)*D
    ratio=lambda z:F((m-p)**2*n,m*(m*n-z))
    removed=ratio(xa)*ratio(xb)*F(m,n)*((m*n-xa-xb)*(1+lower)+F(xa*xb,1)/cf)
    raw=removed+2*h*g*(1+xi)+h*h*(1+chi)
    assert norm==raw
    assert norm>=0,(N,p,m,a,b,norm)
    cross=eta*(1+chi)+g*D
    return norm,cross,1+chi


def exact_gram(N,p,m):
    omega=weights(N,m,p)
    gram=[[F(0) for _ in range(2)] for _ in range(2)]
    for a in range(p+1):
        for b in range(p+1):
            aa,ba,bb=pair(N,p,m,a,b);w=omega[a]*omega[b]
            gram[0][0]+=w*bb;gram[0][1]+=w*ba;gram[1][1]+=w*aa
    gram[1][0]=gram[0][1]
    return gram
