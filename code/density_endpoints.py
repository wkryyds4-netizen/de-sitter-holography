"""Gamma endpoint orders and shifted packets, with every hard fold retained."""
import math
import numpy as np
from scipy.special import gammaincc,digamma
from moving_threshold_curve import ROOT,f,roots,rule,knots,prefactor,fourier,cj,grid

MODELS=[dict(name='gamma1',p=1,beta=1.,gap=0.),
        dict(name='gamma2',p=2,beta=1.,gap=0.),
        dict(name='gap1',p=0,beta=1.,gap=1.)]

def survival(q,p,beta,gap):
    return gammaincc(p+1,beta*np.maximum(0.,np.asarray(q)-gap))

def moments(turn,chi,n,radius,beta,pmax):
    values,restore=np.unique(turn,return_inverse=True)
    edges=knots(radius);x,w=roots(n)
    lo=np.broadcast_to(edges[:-1],(len(values),len(edges)-1)).copy()
    hi=np.maximum(lo,np.minimum(edges[1:][None,:],values[:,None]))
    v=(lo+hi)[:,:,None]/2+(hi-lo)[:,:,None]*x/2
    base=f(v,chi)*np.exp(beta*v)*(hi-lo)[:,:,None]*w/2
    lower=np.array([np.sum(base*v**j,axis=(1,2))[restore] for j in range(pmax+1)])
    hi=np.broadcast_to(edges[1:],lo.shape).copy()
    lo=np.minimum(hi,np.maximum(edges[:-1][None,:],values[:,None]))
    v=(lo+hi)[:,:,None]/2+(hi-lo)[:,:,None]*x/2
    upper=np.sum(f(v,chi)*(hi-lo)[:,:,None]*w/2,axis=(1,2))[restore]
    return lower,upper

def jq_family(chi,n=24,radius=10.,models=MODELS):
    result={}
    groups=sorted(set((m['beta'],m['gap']) for m in models))
    for beta,gap in groups:
        selected=[m for m in models if (m['beta'],m['gap'])==(beta,gap)]
        pmax=max(m['p'] for m in selected)
        u,wu=rule(-radius,radius,n,list(knots(radius))+[-chi/2-gap,chi/2])
        a=np.maximum(0.,-u-chi/2-gap);b=np.maximum(-chi/2,-u)-gap
        low,high=moments(b-a,chi,n,radius,beta,pmax)
        partial_sum=np.zeros(len(u),complex)
        for p in range(pmax+1):
            shifted=sum(math.comb(p,k)*b**(p-k)*(-1)**k*low[k] for k in range(p+1))
            partial_sum+=beta**p/math.factorial(p)*shifted
            inner=np.exp(-beta*b)*partial_sum+gammaincc(p+1,beta*a)*high
            value=np.dot(wu*f(u,chi),inner)
            for m in selected:
                if m['p']==p:result[m['name']]=value
    return result

def jq_halfplane(chi,model,n=24,radius=10.):
    """Independent region restriction and direct survival, no polynomial moments."""
    p,beta,gap=(model[k] for k in ['p','beta','gap'])
    u,wu=rule(-radius,radius,n,list(knots(radius))+[-gap/2,chi/2,-chi/2-gap])
    d=np.where(u>=chi/2,chi/2,u);turn=np.minimum(u,-d-gap)
    edges=knots(radius);x,w=roots(n)
    lo=np.broadcast_to(edges[:-1],(len(u),len(edges)-1)).copy()
    hi=np.maximum(lo,np.minimum(edges[1:][None,:],turn[:,None]))
    v=(lo+hi)[:,:,None]/2+(hi-lo)[:,:,None]*x/2
    lower=np.sum(f(v,chi)*survival(-d[:,None,None]-v,p,beta,gap)*(hi-lo)[:,:,None]*w/2,axis=(1,2))
    lo=np.maximum(edges[:-1][None,:],(-d-gap)[:,None])
    hi=np.maximum(lo,np.minimum(edges[1:][None,:],u[:,None]))
    v=(lo+hi)[:,:,None]/2+(hi-lo)[:,:,None]*x/2
    upper=np.sum(f(v,chi)*(hi-lo)[:,:,None]*w/2,axis=(1,2))
    return 2*np.dot(wu*f(u,chi),lower+upper)

def plateaux(model,n=32,radius=10.):
    p,beta,gap=(model[k] for k in ['p','beta','gap'])
    x,wx=rule(0.,radius,n,list(knots(radius))+[gap])
    kp=f(x,0).real/(16*np.pi**3);km=f(-x,0).real/(16*np.pi**3)
    gain=np.dot(wx,kp);loss=np.dot(wx,km*survival(x,p,beta,gap))
    mix=[];ll=[];diff=[]
    for xx in x:
        u,wu=rule(0.,radius,n,list(knots(radius))+[xx,gap,xx+gap,gap-xx])
        ku=f(-u,0).real/(16*np.pi**3)
        s1=survival(np.maximum(0.,u-xx),p,beta,gap)
        s0=survival(u,p,beta,gap)
        mix.append(np.dot(wu,ku*s1));diff.append(np.dot(wu,ku*(s1-s0)))
        ll.append(np.dot(wu,ku*survival(xx+u,p,beta,gap)))
    mixed=np.dot(wx*kp,mix);lossloss=np.dot(wx*km,ll)
    norm=gain**2+gain*loss+mixed+lossloss
    cross=gain**2+2*gain*loss+lossloss
    edge=2*np.dot(wx*kp,diff)
    return dict(N=float(norm),Xinf=float(cross),Dedge=float(edge),identity_error=float(abs(2*(norm-cross)-edge)))

def coefficients(p,beta,nu=.6):
    F=f(0,0,nu=nu).real
    H=2*digamma(.25+1j*nu/2).real+np.euler_gamma+2*np.log(2)
    d=(np.pi*p-beta*(p+1))/(p+4)
    S=beta**(p+1)*F**2/(512*np.pi**7)
    return dict(F=F,H=float(H),d=float(d),S=float(S))

def tail(model,T):
    p,beta,gap=(model[k] for k in ['p','beta','gap'])
    T=np.asarray(T,dtype=float)
    if gap>0:
        assert p==0
        C=beta*prefactor(-gap)*f(-gap/2,-gap)**2
        return C*np.exp(-1j*gap*T)/T**4
    c=coefficients(p,beta)
    return c['S']*(-1j/(p+3)/(1j*T)**(p+3)+(-c['H']-1j*c['d'])/(1j*T)**(p+4))

def corner_scaled(h,p,beta,gap=0.,n=20):
    """Delta J divided by h^(p+3), evaluated directly on a fixed cube."""
    x,w=roots(n);x=(x+1)/2;w=w/2
    rho,s,eta=np.meshgrid(x,x,x,indexing='ij');wr,ws,we=np.meshgrid(w,w,w,indexing='ij')
    chi=-gap-h;r=gap+h*rho
    a=-r-chi/2;L=r+chi;t=(1-s)*eta
    integrand=(beta**(p+1)/math.factorial(p)*rho**p*np.exp(-beta*h*rho)*(1-rho)**2
               *(1-s)*f(a+L*s,chi)*f(a+L*t,chi))
    return np.sum(wr*ws*we*integrand)
