"""Gaussian coordinate formula for a bounded paired-phase heat-tail diagnostic.

QMC values are empirical estimates of an exact upper-bound functional,
not certified upper bounds and not a test of the pentagon identity.
"""
import numpy as np
from scipy.special import ndtri
from scipy.stats import qmc
from unitary_flip import phase

def normal_nodes(power, seed):
    points=qmc.Sobol(4,scramble=True,seed=seed).random_base2(power)
    return ndtri(points)

def phase_loss(normals, t, h=None, chirp=True):
    u=normals[:,:2]*np.sqrt(np.pi/2)
    d=normals[:,2:]*np.sqrt(2*np.pi**2*t/(1+np.pi*t))
    left=u+d/2;right=u-d/2
    ratio=np.ones(len(u),complex);tails=[]
    if h is not None:
        lp,lt=phase(left[:,0],left[:,1],h)
        rp,rt=phase(right[:,0],right[:,1],h)
        ratio=np.conj(lp)*rp;tails=[lt,rt]
    if chirp:
        ratio*=np.exp(-1j*(u[:,0]*d[:,1]+u[:,1]*d[:,0])/np.pi)
    # |r(s)-r(t)|^2/2 avoids cancellation in 1-Re(ratio).
    loss=abs(1-ratio)**2/2
    return loss,tails

def heat_loss_from_mean(t,mean):
    # Exactly 1-C_t; equality applies to the true expectation, not a QMC estimate.
    return (np.pi*t+mean)/(1+np.pi*t)

def vacuum_heat_error(t):
    a=np.pi*t
    return np.sqrt(2*a*a/((1+a)*(1+2*a)))
