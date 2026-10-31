"""G=0 dS3 group pairing in the declared KAK Haar normalization.

K Haar probability; radial push-forward sinh(r)^2 dr. The linear-gluing
state is a specified model state, not a finite-kappa CLP identification.
"""
import mpmath as mp
from ds3_local_vertex import cosh_denominator,spherical_vertex,vertex
from ds3_mass_spectrum import mu

def spherical(nu,r):
    return mp.mpf(1) if r==0 else mp.sin(nu*r)/(nu*mp.sinh(r))

def pair_schur(nu):return mp.pi/(2*nu**2)
def group_density(nu):return 1/pair_schur(nu)

def triple_tanh(nus):
    S=sum(nus)
    return mp.pi/(8*mp.fprod(nus))*(mp.fsum(mp.tanh(mp.pi*(S-2*n)/2) for n in nus)-mp.tanh(mp.pi*S/2))

def triple_group(nus):
    return mp.pi/(8*cosh_denominator(nus))*mp.fprod(mp.sinh(mp.pi*n)/n for n in nus)

def triple_direct(nus,T=65,step=2):
    grid=list(range(0,T,step))+[T]
    value=mp.quad(lambda r:mp.fprod(mp.sin(n*r) for n in nus)/mp.sinh(r) if r else 0,grid)/mp.fprod(nus)
    tail=2*mp.atanh(mp.exp(-T))/mp.fprod(nus)
    return value,tail

def f_linear(nu):return mp.sech(mp.pi*nu)
def state_density(nu):return group_density(nu)*f_linear(nu)**2
def state_norm(s):
    q=mp.exp(-2*mp.pi*s)
    return 4*s*s/mp.pi**2*q/(1+q)+4*s/mp.pi**3*mp.log1p(q)-2/mp.pi**4*mp.polylog(2,-q)

def emission000(a,b,nu):
    return mp.sqrt(mu(a)*mu(b))*spherical_vertex([a,b,nu],(True,False,True))

def joint(a,b,nu):
    return mu(a)*mu(b)/(8*cosh_denominator([a,b,nu]))

def local_n(a,b):return 1/((1+mp.exp(-2*mp.pi*a))*(1+mp.exp(-2*mp.pi*b)))

def mean_correction(M,w,nu):
    a=M+w/2;b=M-w/2
    A=mp.cosh(2*mp.pi*M);B=mp.cosh(mp.pi*w);C=mp.cosh(mp.pi*nu)
    return (1-w*w/(4*M*M))*(A-B)/(A+C)*(A+1)/(A+B)

def rate(b,nu,s):
    return mp.cosh(mp.pi*b)**2/group_density(b)*mp.quad(
        lambda a:joint(a,b,nu),[s,max(s,b),max(s,b)+2,mp.inf])

def rate_bound(b,s):return b/(2*mp.pi)+1/(2*mp.pi**2)+1/(4*mp.pi**3*s)

def square_tail(L,s):
    beta=2*mp.pi
    return 2/mp.pi**2*mp.exp(-beta*L)*((L**3-s*s*L)/beta+
        (3*L*L-s*s)/beta**2+6*L/beta**3+6/beta**4)

def compact_sine(r,a,b):
    # Exact endpoint integration of (nu-a)^4 (b-nu)^4 exp(i nu r).
    L=b-a
    if abs(r)<2:return mp.quad(lambda nu:(nu-a)**4*(b-nu)**4*mp.sin(nu*r),[a,b])
    coeff=[mp.mpf(0)]*9
    for n in range(5):coeff[n+4]=(-1)**n*mp.binomial(4,n)*L**(4-n)
    terms=[]
    for j in range(9):
        end=mp.fsum(coeff[n]*mp.factorial(n)/mp.factorial(n-j)*L**(n-j) for n in range(j,9))
        start=coeff[j]*mp.factorial(j)
        terms.append((-1)**j*(end*mp.exp(1j*b*r)-start*mp.exp(1j*a*r))/(1j*r)**(j+1))
    return mp.im(mp.fsum(terms))

def compact_pair_norm(a,b):return mp.pi/2*(b-a)**17*mp.beta(9,9)
def compact_pair_tail(R,a,b):
    coeff=[24,-480,2160,-3360,1680]
    integral=mp.fsum(mp.mpf(x*y)/(i+j+1) for i,x in enumerate(coeff) for j,y in enumerate(coeff))
    return mp.pi/2*(b-a)**9*integral/R**8
