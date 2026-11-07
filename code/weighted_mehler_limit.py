"""Pure exact-product/modular forms for fixed integer Mehler channels."""
import mpmath as mp

def leading_gap(z,d):
    base=mp.mpf(1)/mp.pi if not z else z/mp.sinh(mp.pi*z)
    return mp.mpf(2)**(2*d-1)/mp.factorial(2*d-1)*base*mp.fprod(j*j+z*z for j in range(1,d))

def leading_density(y,z,d):
    return mp.exp(-y*y/2)/mp.sqrt(2*mp.pi)*leading_gap(z,d)

def image_quotient(x,h,terms=8):
    """T(x)/[Gaussian*(1-exp(-2 pi b/h))], with endpoint limits."""
    b=min(x,2*mp.pi-x)
    assert b>=0
    a=mp.exp(-2*mp.pi*b/h)
    return mp.fsum([mp.mpf(1)]+[
        (-1)**n*mp.exp(-2*mp.pi*n*((n+1)*mp.pi-b)/h)
        *mp.fsum(a**j for j in range(2*n+1)) for n in range(1,terms+1)])

def image_tail(h,terms=8):
    n=terms+1;c=2*mp.pi**2/h;r=mp.exp(-c*(2*n+1))
    return mp.exp(-c*n*n)*((2*n+1)/(1-r)+2*r/(1-r)**2)

def modular_density(h,y,z,d=2,terms=8):
    h,y,z=map(mp.mpf,(h,y,z));z=abs(z)
    a=mp.sqrt(h)*y;delta=h*z;S=mp.pi+a
    theta=(S+delta)/2;phi=(S-delta)/2
    assert 0<theta<mp.pi and 0<phi<mp.pi
    q=mp.exp(-h);dual=mp.exp(-4*mp.pi**2/h)
    D=mp.fprod(-mp.expm1(-j*h) for j in range(1,2*d))
    def edge(x):
        b=min(x,2*mp.pi-x)
        return -mp.expm1(-2*mp.pi*b/h)
    images=(edge(2*theta)*edge(2*phi)/edge(S)
        *image_quotient(2*theta,h,terms)*image_quotient(2*phi,h,terms)
        /image_quotient(S,h,terms)/image_quotient(delta,h,terms))
    ratio=h/(2*mp.pi) if not z else mp.sin(delta/2)/mp.sinh(mp.pi*z)
    factors=mp.fprod(4*q**(2*j)*(mp.cosh(j*h)+mp.cos(a))
                     *(mp.cosh(j*h)-mp.cos(delta)) for j in range(1,d))
    return ((2*mp.pi)**mp.mpf('1.5')/mp.pi**2/D
        *mp.exp(h/8-y*y/2-h*z*z/2)*mp.qp(dual,dual)**3
        *mp.sin(theta)*mp.sin(phi)*mp.sin(S/2)*ratio*factors*images)

def product_density(h,y,z,d=2):
    h,y,z=map(mp.mpf,(h,y,z));q=mp.exp(-h);r=q**d
    theta=(mp.pi+mp.sqrt(h)*y+h*z)/2
    phi=(mp.pi+mp.sqrt(h)*y-h*z)/2
    def A(x):return abs(mp.qp(mp.exp(1j*x),q))**2
    denominator=abs(mp.qp(r*mp.exp(1j*(theta+phi)),q)
                    *mp.qp(r*mp.exp(1j*(theta-phi)),q))**2
    return h**mp.mpf('1.5')/2*mp.qp(q,q)**2/(4*mp.pi**2)*A(2*theta)*A(2*phi)*mp.qp(r*r,q)/denominator

def coefficient(y,z,d=2):
    return mp.mpf(4*d+1)/8-mp.mpf(2*d+1)/8*y*y-z*z/2

def symmetric_weight(h,y,z):
    return 2*(1-mp.cos(mp.sqrt(h)*y)*mp.cos(h*z))/(-mp.expm1(-h))-1

def gap(h,y,z):
    return -4*mp.cos(mp.sqrt(h)*y/2)*mp.sin(h*z/2)/mp.sqrt(h*(-mp.expm1(-h)))

def response(u,d=2):
    return mp.mpf(d)/2/mp.cosh(u)**(2*d)*(u*mp.tanh(u)-(1+2*d)*mp.tanh(u)**2)
