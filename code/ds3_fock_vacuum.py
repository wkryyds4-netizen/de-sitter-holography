"""Specified adjoint-completed half-period field: vacuum return block.

No finite-G or unsmeared-field claim. Full CG/Fock existence is analytic,
not established by the finite-dimensional auxiliary diagnostics.
"""
import mpmath as mp
from ds3_two_matter import tilted,tilt_constant,lr
from ds3_group_pairing import emission000,group_density,triple_group

def energy(x):return mp.sqrt(1+x*x)

def block_bound(s):
    s=mp.mpf(s)
    return tilt_constant(s)*(1/(2*mp.pi**2)+1/(2*s*mp.pi**3))

def raw_rate(a,b,p):
    E=emission000(a,b,p)
    w=group_density(b)*triple_group([a,p,b])
    return abs(E)**2*mp.exp(mp.pi*(energy(a)-energy(b)))/w

def one_step_tail(D,s):
    D,s=map(mp.mpf,[D,s]);P=mp.pi
    return tilt_constant(s)*mp.exp(-P*D)/(2*P)*(1/P+(D/P+1/P**2)/s)

def return_block(M,p,s,t,D=24):
    M,p,s,t,D=map(mp.mpf,[M,p,s,t,D])
    if M<s:raise ValueError('Initial label below threshold')
    lo=max(s-M,-D)
    grid=sorted(set([lo,D]+[mp.mpf(j) for j in range(-int(D),int(D)+1,2) if lo<j<D]+([mp.mpf(0)] if lo<0<D else [])))
    def integrand(u):
        a=M+u
        delta=u*(a+M)/(energy(a)+energy(M))
        return tilted(a,M,p)*mp.exp(-1j*t*delta)
    value=mp.quad(integrand,grid)
    limit=lr(t,p);first=1j*mp.diff(lambda tau:lr(tau,p),t)/M
    return dict(M=M,p=p,s=s,t=t,D=D,value=value,tail=one_step_tail(D,s),
        limit=limit,first_correction=first,second_scaled=M*M*(value-limit-first),
        vacuum_commutator_norm=4*mp.im(value)**2,
        predicted_scaled_vacuum=4*mp.diff(lambda tau:lr(tau,p),t)**2)
