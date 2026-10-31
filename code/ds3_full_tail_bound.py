"""Exact rational ratio certificates and exact-arithmetic full scalar-kernel tails.

Finite sums and Gamma coefficients use mpmath, NOT interval arithmetic.
The polynomial certificate is exact; the evaluated bounds exclude rounding.
No heavy-mass uniformity or distributional convergence is claimed.
"""
from functools import lru_cache
import sympy as sp
import mpmath as mp
from global_racah import shadow_coefficient
from ds3_crossing_dictionary import period


def numeric(z):
    re, im = sp.expand(z).as_real_imag()
    def rat(x):
        n, d = x.as_numer_denom()
        return mp.mpf(int(n))/int(d)
    return mp.mpc(rat(re), rat(im))


@lru_cache(maxsize=None)
def certificate(top, bottom):
    """Prove |t[n+1]/t[n]| <= ((n+1)/(n+2))**(3/2), n>=N0.

    Rational imaginary parts only. All shifted polynomial coefficients
    are checked exactly, not just the polynomial at a sampling grid.
    """
    assert len(top) == 4 and len(bottom) == 3
    assert all(sp.re(a) == sp.Rational(1, 2) for a in top)
    assert all(sp.re(b) == 1 for b in bottom)
    x = sp.Symbol('x', real=True)
    numerator = sp.prod((x+sp.Rational(1, 2))**2+sp.im(a)**2 for a in top)
    denominator = (x+1)**2*sp.prod((x+1)**2+sp.im(b)**2 for b in bottom)
    poly = sp.Poly(denominator*(x+1)**3-numerator*(x+2)**3, x)
    assert poly.degree() == 10 and poly.LC() == 1
    N = 1
    while True:
        shifted = poly.shift(N)
        coeffs = shifted.all_coeffs()
        if all(c >= 0 for c in coeffs):
            return dict(N0=N, degree=poly.degree(),
                        coefficients=[str(c) for c in poly.all_coeffs()],
                        shifted_coefficients=[str(c) for c in coeffs])
        N *= 2
        if N > 2**24:
            raise ValueError('Certificate search exceeded declared cost cap')


def series(top, bottom, N, weighted=False):
    cert = certificate(tuple(top), tuple(bottom))
    if N < cert['N0']:
        raise ValueError('Cutoff lies before exact polynomial certificate')
    a, b = list(map(numeric, top)), list(map(numeric, bottom))
    b4 = b+[mp.mpc(1)]
    t = mp.mpc(1)
    q = sum(mp.digamma(v) for v in b4)-sum(mp.digamma(v) for v in a) if weighted else None
    values = []
    for n in range(N):
        values.append(t*q if weighted else t)
        t *= mp.fprod(n+v for v in a)/(mp.fprod(n+v for v in b)*(n+1))
        if weighted:
            q += sum(1/(n+v) for v in b4)-sum(1/(n+v) for v in a)
    if weighted:
        L = sum(abs(v-w) for v, w in zip(b4, a))
        c = L*(1+1/(mp.mpf(N)+mp.mpf('.5')))*(N+1)/(mp.mpf(N)+mp.mpf('.5'))
        tail = abs(t)*c*(1/mp.mpf(N+1)+mp.mpf(2)/3)
    else:
        L = None
        tail = abs(t)*(2*N+3)
    value = mp.fsum(values)
    return value, tail, dict(N=N, first_omitted=t, weighted=weighted, L=L,
                            certificate=cert, absolute_sum=mp.fsum(abs(v) for v in values))


def parameters(H, hp, hs):
    h1, h2, h3, h4 = hs
    hp0 = sp.expand(H+h1+h3-1)
    z = sp.expand(hp0-hp)
    top = tuple(map(sp.expand, [H+h1+h2-1,H+h3+h4-1,H+h3-h4,H+h1-h2]))
    if z == 0:
        return z, [(top, tuple(map(sp.expand,[2*hp0,2*H,sp.Integer(1)])))]
    return z, [(tuple(sp.expand(a-z) for a in top),tuple(map(sp.expand,[2*hp,2*H-z,1-z]))),
               (top,tuple(map(sp.expand,[2*hp0-z,2*H,1+z])))]


def omega(H, hp, hs, refinement=1):
    z, groups = parameters(H, hp, hs)
    certs = [certificate(a, b) for a, b in groups]
    N = refinement*max(64, *[c['N0'] for c in certs])
    if N > 65536:
        raise ValueError('Direct-sum cost cap exceeded')
    pieces = [series(a,b,N,weighted=(z==0)) for a,b in groups]
    if z == 0:
        value, tail, record = pieces[0]
        return value, tail, dict(collision=True, groups=[record])
    h1,h2,h3,h4 = hs
    hp0 = sp.expand(H+h1+h3-1)
    A = list(map(numeric,[H+h1+h2-1,H+h3+h4-1,H+h3-h4,H+h1-h2]))
    H,hp,hp0,z = map(numeric,[H,hp,hp0,z])
    ca = mp.gamma(z)*mp.gamma(2*H)*mp.gamma(A[3]-z)*mp.gamma(A[2]-z)/(mp.gamma(A[3])*mp.gamma(A[2])*mp.gamma(2*H-z))
    cb = mp.gamma(-z)*mp.gamma(2*hp)*mp.gamma(A[0])*mp.gamma(A[1])/(mp.gamma(A[0]-z)*mp.gamma(A[1]-z)*mp.gamma(2*hp0-z))
    value = ca*pieces[0][0]+cb*pieces[1][0]
    tail = abs(ca)*pieces[0][1]+abs(cb)*pieces[1][1]
    return value, tail, dict(collision=False, z=z, coefficients=[ca,cb],
                            groups=[p[2] for p in pieces],
                            pair_condition=(abs(ca*pieces[0][0])+abs(cb*pieces[1][0]))/max(abs(value),mp.mpf('1e-1000')))


def physical(b,d,external,refinement=1):
    b,d = map(sp.Rational,[str(b),str(d)])
    external = [sp.Rational(str(v)) for v in external]
    hs = [(1+sp.I*v)/2 for v in external]
    H=(1+sp.I*b)/2;hp=(1+sp.I*d)/2
    hn=list(map(numeric,hs));h1,h2,h3,h4=hn;h=numeric(H)
    pref=mp.gamma(h+h1-h2)*mp.gamma(h-h1+h2)*mp.gamma(1-h+h4-h3)*mp.gamma(1-h+h3-h4)/(4*mp.gamma(2*h)*mp.gamma(2-2*h))
    terms=[];errors=[];records=[]
    for hq,co in [(hp,shadow_coefficient(numeric(1-hp),h1,h4)),(1-hp,shadow_coefficient(numeric(hp),h3,h2))]:
        parts=[omega(x,hq,hs,refinement) for x in [H,1-H]]
        q=numeric(hq)
        k=co*pref*mp.sin(mp.pi*(q-h1-h4))*mp.sin(mp.pi*(q-h2-h3))
        x,ex,_=parts[0];y,ey,_=parts[1]
        terms.append(k*x*y)
        errors.append(abs(k)*(abs(x)*ey+abs(y)*ex+ex*ey))
        records.append(dict(coefficient=k,omega=[dict(value=v,tail=e,details=r) for v,e,r in parts]))
    a,q,c,p=map(numeric,external);bn,dn=numeric(b),numeric(d)
    ns=period(a,q,bn)*period(c,p,-bn);nt=period(c,q,dn)*period(a,p,-dn)
    factor=bn*dn/mp.pi**4*(ns/abs(ns))/(nt/abs(nt))
    return factor*mp.fsum(terms),abs(factor)*mp.fsum(errors),dict(shadows=records,normalization=factor,
            shadow_condition=mp.fsum(abs(v) for v in terms)/max(abs(mp.fsum(terms)),mp.mpf('1e-1000')))
