"""Exact four-leg Gaussian filter covariance and conditional eikonal comparison."""
import numpy as np
from scipy.special import ndtri
from scipy.stats import qmc
from principal_eikonal import kernel

def covariance(delta):
    return delta**2*np.array([[5,-1,-3,-1],[-1,5,-1,-3],[-3,-1,5,-1],[-1,-3,-1,5]])/8

def slope(n,nu,nup,center_variance_ratio=.5):
    r=lambda k,v:((k+2)**2/4+v*v)/(k+3)
    return center_variance_ratio*n*n/2-r(n,nu)+r(0,nu)-r(n,nup)+r(0,nup)

def paired_integral(delta,u,power,seed,base,m=1000,gamma=.25,nu=.6,nup=.9):
    sampler=qmc.Sobol(8,scramble=True,seed=seed);sums=np.zeros(6,complex)
    for _ in range(2**(power-15)):
        a=np.clip(sampler.random(2**15),np.finfo(float).eps,1-np.finfo(float).eps)
        normals=ndtri(a[:,:3]);X=delta*normals[:,0]/np.sqrt(2);Y=delta*normals[:,1]/np.sqrt(2)
        r=-np.log(a[:,3:7]);A=np.exp(Y)*np.cosh(r[:,0])+np.exp(-Y)*np.cosh(r[:,1])
        B=np.exp(X)*np.cosh(r[:,2])+np.exp(-X)*np.cosh(r[:,3])
        weight=np.exp(r.sum(axis=1))/(A*A*B*B)
        weight*=np.cos(nu*r[:,0])*np.cos(nu*r[:,1])*np.cos(nup*r[:,2])*np.cos(nup*r[:,3])
        phi=np.pi*a[:,7];h=-np.cos((1-gamma)*(np.pi-phi))/(2*m*np.sinc(gamma))
        values=[]
        for scale in [1/np.sqrt(2),1.]:
            T=np.log(m)+u+delta*scale*normals[:,2]
            values.append(weight*(kernel(h*np.exp(T)/(A*B))+kernel(h*np.exp(-T)/(A*B),-1)-2))
        sums+=np.array([v.sum() for v in values]+[(values[0]-values[1]).sum(),weight.sum(),
                       np.abs(weight+values[0]).sum(),np.abs(values[0]-values[1]).sum()])
    means=sums/(2**power*base)
    pair=lambda v:[float(v.real),float(v.imag)]
    return dict(delta=delta,u=u,power=power,seed=seed,m=m,gamma=gamma,nu=nu,nup=nup,base=base,
                correlated=pair(1+means[0]),independent=pair(1+means[1]),difference=pair(means[2]),
                baseline_ratio=float(means[3].real),raw_correlated=pair(means[0]+means[3]),
                absolute_correlated=float(means[4].real),absolute_difference=float(means[5].real))
