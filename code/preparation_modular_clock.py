"""Pure spectral tools for the selected state's modular clock, at fixed q."""
import numpy as np
from scipy.linalg import eigh_tridiagonal
from numpy.polynomial.legendre import leggauss

def parameters(q):
    lam = -np.log(q)/4
    b = np.sqrt(q)*(1+lam)
    a = 1-b
    R = 2/np.sqrt(1-q)
    return lam, a, b, R

def angular_density(theta, q, products):
    logmu = np.log(2/np.pi)+2*np.log(np.sin(theta))
    for n in range(1, products+1):
        z = q**n
        logmu += np.log1p(-z)+np.log(1-2*z*np.cos(2*theta)+z*z)
    return np.exp(logmu)

def mehler(theta, phi, q, r, products):
    cp = np.cos(theta[:, None]+phi[None, :])
    cm = np.cos(theta[:, None]-phi[None, :])
    lk = np.zeros_like(cp)
    for n in range(products):
        z = r*q**n
        lk += np.log1p(-r*r*q**n)-np.log(1-2*z*cp+z*z)-np.log(1-2*z*cm+z*z)
    return np.exp(lk)

def angular(q, L, products=192):
    z, w = leggauss(L)
    theta = (z+1)*np.pi/2
    mu = w*np.pi/2*angular_density(theta, q, products)
    E = 2*np.cos(theta)/np.sqrt(1-q)
    return E, mu, mu[:, None]*mu[None, :]*mehler(theta, theta, q, q, products)

def jacobi(q, L):
    E, V = eigh_tridiagonal(np.zeros(L), np.sqrt((1-q**np.arange(1,L))/(1-q)))
    V *= np.where(V[0]<0, -1., 1.)[None, :]
    root = V[0]
    transfer = V.T@(q**np.arange(L)[:, None]*V)
    return E, root*root, root[:, None]*root[None, :]*transfer

def spectrum(q, data):
    E, mu, joint = data
    lam, a, b, R = parameters(q)
    rho = a+b*E*E
    k = -np.log(rho)
    delta = k[:, None]-k[None, :]
    weights = rho[:, None]*joint
    return dict(E=E, mu=mu, joint=joint, rho=rho, delta=delta, weights=weights,
                lam=lam, a=a, b=b, R=R)

def correlation(spec, s):
    terms = spec['weights']*np.exp(1j*s*spec['delta'])
    value = terms.sum()
    return value, float(np.abs(terms).sum()/max(abs(value), 1e-300))

def sum_rules(spec):
    weights, delta = spec['weights'], spec['delta']
    return dict(mass=float(weights.sum()),
                first_frequency=float(np.sum(weights*(-delta))),
                second_frequency=float(np.sum(weights*delta**2)),
                relative_entropy=float(np.sum(spec['mu']*spec['rho']*np.log(spec['rho']))),
                min_joint_weight=float(spec['joint'].min()))

def product_log_error(q, r, products):
    # Uniform absolute log error for two mu factors and one Mehler kernel.
    return (6*q**(products+1)/(1-q**(products+1))+
            r*r*q**products/(1-r*r*q**products)+
            4*r*q**products/(1-r*q**products))/(1-q)

def window_rates(q, output_interval, input_interval, L=64, products=192):
    _, _, _, R = parameters(q)
    z, w = leggauss(L)
    def grid(interval):
        lo, hi = np.arccos(interval[1]/R), np.arccos(interval[0]/R)
        theta = lo+(z+1)*(hi-lo)/2
        E = R*np.cos(theta)
        mu = w*(hi-lo)/2*angular_density(theta, q, products)
        return theta, E, mu
    th, E, mx = grid(output_interval)
    ph, F, my = grid(input_interval)
    mass = mx[:, None]*my[None, :]*mehler(th, ph, q, q*q, products)
    f = lambda x: 1-np.sqrt(q)+np.sqrt(q)*x*x
    forward = np.sum(mass*f(F)[None, :])
    reverse = np.sum(mass*f(E)[:, None])
    return dict(forward=float(forward), reverse=float(reverse),
                ratio=float(reverse/forward), total_joint_mass=float(mass.sum()))
