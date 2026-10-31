"""Real Gaussian smearing of the original local field, distinct from balancing."""
import mpmath as mp
from ds3_two_matter import tilted,k
from ds3_fock_vacuum import energy,one_step_tail,block_bound

def multiplier_bound(sigma):return mp.exp(mp.pi**2/(4*mp.mpf(sigma)**2))

def density(u,p,sigma):return mp.exp(-sigma*sigma*u*u-mp.pi*u)*k(u,p)

def limiting(t,p,sigma,D=32,derivative=0):
    t,p,sigma=mp.mpc(t),mp.mpf(p),mp.mpf(sigma)
    return mp.quad(lambda u:(-1j*u)**derivative*density(u,p,sigma)*mp.exp(-1j*t*u),list(range(-D,D+1,2)))

def return_block(M,p,s,t,sigma,D=32):
    M,p,s,t,sigma,D=map(mp.mpf,[M,p,s,t,sigma,D]);lo=max(s-M,-D)
    grid=sorted(set([lo,D]+[mp.mpf(j) for j in range(-int(D),int(D)+1,2) if lo<j<D]+([mp.mpf(0)] if lo<0<D else [])))
    def fun(u):
        a=M+u;delta=u*(a+M)/(energy(a)+energy(M))
        return tilted(a,M,p)*mp.exp(-sigma*sigma*delta*delta-mp.pi*delta-1j*t*delta)
    value=mp.quad(fun,grid)
    return dict(M=M,p=p,s=s,t=t,sigma=sigma,D=D,value=value,
        tail=multiplier_bound(sigma)*one_step_tail(D,s),vacuum_commutator_norm=4*mp.im(value)**2)
