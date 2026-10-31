"""Wilson-duality reflection of LPRS Omega and a two-factor physical kernel.

Uses a meromorphic special-function identity, not Wilson-transform unitarity.
The physical complex phases are retained; no empirical real projection.
"""
import sympy as sp
import mpmath as mp
from ds3_full_tail_bound import numeric, omega
from global_racah import shadow_coefficient
from ds3_crossing_dictionary import period


def log_reflection(H,hp,hs):
    h1,h2,h3,h4=hs
    top=[H+h1+h2-1,H+h3+h4-1,hp+1-h1-h4,hp+1-h2-h3]
    bottom=[H+1-h1-h2,H+1-h3-h4,hp+h1+h4-1,hp+h2+h3-1]
    return sum(mp.loggamma(v) for v in top)-sum(mp.loggamma(v) for v in bottom)


def reflection(H,hp,hs):
    return mp.exp(log_reflection(H,hp,hs))


def wilson_data(H,hp,hs):
    h1,h2,h3,h4=hs
    hp0=H+h1+h3-1;z=hp0-hp
    A=[H+h1+h2-1,H+h3+h4-1,H+h3-h4,H+h1-h2]
    T=[v-z/2 for v in A]
    lam=z/(2j);x=(hp0-H-z/2)/1j
    K=mp.fprod(mp.gamma(v) for v in [2*H,2*hp,A[0],A[1],A[2]-z,A[3]-z])
    return T,lam,x,K


def physical(b,d,external,refinement=1):
    b,d=map(sp.Rational,[str(b),str(d)])
    external=[sp.Rational(str(v)) for v in external]
    hs=[(1+sp.I*v)/2 for v in external];H=(1+sp.I*b)/2;hp=(1+sp.I*d)/2
    parts=[omega(x,hp,hs,refinement) for x in [H,1-H]]
    hn=list(map(numeric,hs));h1,h2,h3,h4=hn;h=numeric(H);q=numeric(hp)
    pref=mp.gamma(h+h1-h2)*mp.gamma(h-h1+h2)*mp.gamma(1-h+h4-h3)*mp.gamma(1-h+h3-h4)/(4*mp.gamma(2*h)*mp.gamma(2-2*h))
    plus=shadow_coefficient(1-q,h1,h4)*pref*mp.sin(mp.pi*(q-h1-h4))*mp.sin(mp.pi*(q-h2-h3))
    minus=shadow_coefficient(q,h3,h2)*pref*mp.sin(mp.pi*(1-q-h1-h4))*mp.sin(mp.pi*(1-q-h2-h3))
    logR=log_reflection(h,1-q,hn)+log_reflection(1-h,1-q,hn)
    # Exponentiation combines the enormous coefficient and tiny reflection product.
    effective=mp.exp(mp.log(minus)+logR)
    a,qm,c,p=map(numeric,external);bn,dn=numeric(b),numeric(d)
    ns=period(a,qm,bn)*period(c,p,-bn);nt=period(c,qm,dn)*period(a,p,-dn)
    factor=bn*dn/mp.pi**4*(ns/abs(ns))/(nt/abs(nt))
    plus*=factor;effective*=factor
    x,ex,_=parts[0];y,ey,_=parts[1];P=x*y
    errorP=abs(x)*ey+abs(y)*ex+ex*ey
    terms=[plus*P,effective*mp.conj(P)]
    value=mp.fsum(terms);tail=(abs(plus)+abs(effective))*errorP
    return value,tail,dict(coefficients=[plus,effective],log_reflection_product=logR,
               original_minus_modulus=abs(factor*minus),coefficient_modulus_residual=abs(abs(effective)/abs(plus)-1),
               omega=[dict(value=v,tail=e,details=r) for v,e,r in parts],
               summation_condition=mp.fsum(abs(t) for t in terms)/max(abs(value),mp.mpf('1e-1000')))
