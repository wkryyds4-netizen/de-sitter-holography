"""de Sitter_4 conformal scalar worldline spectra with explicit modular units.

W(t)=int rho(E) exp(-i E s) dE, s=t/beta, beta=2*pi*L.
All reported correlators use A=phi_f/sqrt(<phi_f^2>).
No gravitational interactions or observer recoil are included.
"""
from dataclasses import dataclass
import math
import numpy as np
from scipy.integrate import quad, quad_vec
from scipy.special import roots_legendre
from clock_model import CompactClock, SmoothClock

def rho_shape(energy):
    E = np.asarray(energy, dtype=float)
    ans = np.empty_like(E)
    small = abs(E) < 1e-5
    positive = (E >= 0) & ~small
    negative = (E < 0) & ~small
    ans[small] = 1 + E[small]/2 + E[small]**2/12 - E[small]**4/720
    ans[positive] = E[positive]/(-np.expm1(-E[positive]))
    x = -E[negative]
    ans[negative] = x*np.exp(-x)/(-np.expm1(-x))
    return ans

def physical_rho(energy, radius=1):
    beta = 2*np.pi*radius
    return rho_shape(energy)/(4*np.pi**2*beta**2)

@dataclass(frozen=True)
class Smearing:
    family: str
    sigma: float

    def __post_init__(self):
        if self.family not in ("gaussian", "cosine") or self.sigma <= 0:
            raise ValueError("family gaussian/cosine and sigma>0 required")

    @property
    def half_duration(self):
        return self.sigma / math.sqrt(1/3 - 2/np.pi**2) if self.family == "cosine" else math.inf

    def fourier(self, E):
        E = np.asarray(E, dtype=float)
        if self.family == "gaussian":
            return np.exp(-.5*self.sigma**2*E**2)
        z = E*self.half_duration/np.pi
        return np.sinc(z)+.5*(np.sinc(z-1)+np.sinc(z+1))

    def density(self, E):
        return rho_shape(E)*self.fourier(E)**2

    def time_function(self, u):
        u = np.asarray(u, dtype=float)
        if self.family == "gaussian":
            return np.exp(-u*u/(2*self.sigma**2))/(math.sqrt(2*np.pi)*self.sigma)
        T = self.half_duration
        return np.where(abs(u)<=T, np.cos(np.pi*u/(2*T))**2/T, 0)

    def tail_bound_unnormalized(self, radius):
        """Rigorous bound on int_{|E|>R} rho_shape |fhat|^2, R positive."""
        R = float(radius)
        if R <= 0:
            return math.inf
        thermal = (1+math.exp(-R))/(-math.expm1(-R))
        if self.family == "gaussian":
            return thermal*math.exp(-(self.sigma*R)**2)/(2*self.sigma**2)
        T = self.half_duration
        if R <= np.pi/T:
            return math.inf
        coefficient = np.pi**4/(T**6*(1-(np.pi/(T*R))**2)**2)
        return thermal*coefficient/(4*R**4)

def composite_nodes(edges, order=16):
    nodes, weights = roots_legendre(order)
    e = np.asarray(edges)
    centers, halves = (e[1:]+e[:-1])/2, (e[1:]-e[:-1])/2
    return (centers[:,None]+halves[:,None]*nodes).ravel(), (halves[:,None]*weights).ravel()

@dataclass
class Spectrum:
    smearing: Smearing
    order: int = 16
    periods: int = 40

    def __post_init__(self):
        if self.smearing.family == "gaussian":
            self.radius = max(36., 8/self.smearing.sigma)
            step = min(2., 1/self.smearing.sigma)
        else:
            step = np.pi/self.smearing.half_duration
            self.radius = max(36., self.periods*step)
        edges = np.unique(np.r_[-self.radius, np.arange(-self.radius,0,step), 0,
                               np.arange(step,self.radius,step), self.radius])
        self.edges = edges
        self.energy, base = composite_nodes(edges, self.order)
        unnorm = base*self.smearing.density(self.energy)
        self.mass_truncated = float(unnorm.sum())
        self.weights = unnorm/self.mass_truncated
        self.tail_bound = self.smearing.tail_bound_unnormalized(self.radius)
        # Numerical mass must be convergence-checked independently.
        self.discarded_fraction_bound = self.tail_bound/self.mass_truncated
        self.physical_variance = self.mass_truncated/(16*np.pi**4)
        self.x = self.energy[:,None]
        self.y = self.energy[None,:]
        self.joint = self.weights[:,None]*self.weights[None,:]
        self._cdf_nodes, self._cdf_weights = roots_legendre(max(32,self.order))
        self._positive_edges = edges[edges>=0]
        e = self._positive_edges
        mid, half = (e[1:]+e[:-1])/2, (e[1:]-e[:-1])/2
        samples = mid[:,None]+half[:,None]*self._cdf_nodes
        masses = half*(self.smearing.density(samples)@self._cdf_weights)/self.mass_truncated
        self._right_mass = np.r_[np.cumsum(masses[::-1])[::-1],0]
        self._tails_cache = {}

    def c2(self, times):
        return np.exp(-1j*np.asarray(times)[...,None]*self.energy)@self.weights

    def q_arrays(self):
        return (np.maximum(np.maximum(self.x,self.y),0),
                np.maximum(self.x,0)+np.maximum(self.y,0),
                np.maximum(0,self.x+np.maximum(self.y,0)))

    def _positive_tail(self, q):
        """U(q)=int_q^R mu for q>=0, integrated as a tail, not 1-CDF."""
        values=np.asarray(q,dtype=float)
        flat=values.ravel()
        out=np.zeros_like(flat)
        active=(flat<self.radius)
        z=np.maximum(0,flat[active])
        idx=np.searchsorted(self._positive_edges,z,side="right")-1
        idx=np.clip(idx,0,len(self._positive_edges)-2)
        right=self._positive_edges[idx+1]
        half=(right-z)/2
        nodes=(right+z)[:,None]/2+half[:,None]*self._cdf_nodes
        partial=half*(self.smearing.density(nodes)@self._cdf_weights)/self.mass_truncated
        out[active]=partial+self._right_mass[idx+1]
        return out.reshape(values.shape)

    def _tail_convolution(self, p, low, high):
        if high<=low:
            return 0.
        step=min(2.,np.pi/self.smearing.half_duration) if self.smearing.family=="cosine" else min(2.,1/self.smearing.sigma)
        edges=np.linspace(low,high,max(1,math.ceil((high-low)/step))+1)
        x,w=composite_nodes(edges,max(24,self.order))
        return float(np.dot(w,self.smearing.density(x)*self._positive_tail(p-x))/self.mass_truncated)

    def _spectral_tails(self, p):
        """Stable CDF tails for the four spectral tasks; no indicator grid."""
        p=float(p)
        if p<0:
            raise ValueError("A positive intermediate-energy threshold is required")
        if p>=2*self.radius:
            return np.zeros(4)
        if p in self._tails_cache:
            return self._tails_cache[p]
        U=float(self._positive_tail(p))
        F0=1-float(self._positive_tail(0))
        convolution=self._tail_convolution(p,max(0,p-self.radius),min(p,self.radius))
        negative=self._tail_convolution(p,max(-self.radius,p-self.radius),0)
        answer=np.array([U,2*U-U*U,(1+F0)*U+convolution,U+convolution+negative])
        self._tails_cache[p]=answer
        return answer

    def tails(self, Q):
        values=self._spectral_tails(Q)
        return {"T2_normalized":float(values[0]),"T4_normalized":float(np.mean(values[1:])),
                "T4_pieces":values[1:].tolist(),"method":"CDF tail convolution, finite spectral radius"}

    def errors_many(self, clocks):
        """Average task tails over clock energy; independent of 2D q-grid errors.

        Gaussian/finite-duration spectral truncation and clock integration tails
        are separate from quadrature error. Inner GL convergence is checked by
        validate_field.py against independent adaptive q-density integrals.
        """
        clocks=list(clocks)
        uppers=[]; clock_bounds=[]; points=[]
        for clock in clocks:
            if hasattr(clock,"maximum_energy"):
                upper=clock.maximum_energy
                points.extend([clock.threshold,upper])
            elif hasattr(clock,"scale"):
                upper=clock.threshold+80*clock.scale
                points.append(clock.threshold)
            else:
                upper=40*clock.mean
            upper=min(upper,2*self.radius)
            uppers.append(upper)
            clock_bounds.append(0. if upper>=2*self.radius else float(clock.kernel(upper)))
        upper=max(uppers)
        step=min(2.,np.pi/self.smearing.half_duration) if self.smearing.family=="cosine" else 2.
        points=np.unique(np.r_[points,uppers,np.arange(0,upper,step),self.radius])
        points=points[(points>0)&(points<upper)]
        def integrand(p):
            densities=np.array([float(c.density(p)) if p<=hi else 0. for c,hi in zip(clocks,uppers)])
            return self._spectral_tails(p)[:,None]*densities[None,:]
        values,estimate=quad_vec(integrand,0,upper,points=points,epsabs=2e-10,epsrel=2e-10,limit=10000)
        records=[]
        for j,clock in enumerate(clocks):
            d2=float(values[0,j]); deltas=values[1:,j]
            c2zero=1-d2; c4zero=3-float(deltas.sum())
            records.append({"sup_C2_absolute":d2,"sup_C4_absolute_over_3":float(np.mean(deltas)),
                "C4_pairing_deficits":deltas.tolist(),"C2_at_zero":c2zero,
                "C4_at_zero":c4zero,"dressed_fourth_cumulant":c4zero-3*c2zero**2,
                "spectral_truncation_C2_bound":2*self.discarded_fraction_bound,
                "spectral_truncation_C4_over_3_bound":4*self.discarded_fraction_bound,
                "clock_integration_tail_bound":clock_bounds[j],
                "outer_quadrature_error_estimate":float(estimate),
                "method":"clock-density average of stable CDF task tails; GL inner convolution"})
        return records

    def errors(self, clock):
        return self.errors_many([clock])[0]

    def c4(self, times, clock=None):
        """Literal ordered product A(t1) A(t2) A(t3) A(t4). No i0 regulator after smearing."""
        t1,t2,t3,t4 = map(float,times)
        phase_pairs = ((t2-t1,t4-t3),(t3-t1,t4-t2),(t4-t1,t3-t2))
        value = 0j
        for q,(a,b) in zip(self.q_arrays(),phase_pairs):
            kernels = 1 if clock is None else clock.kernel(q)
            value += np.sum(self.joint*kernels*np.exp(1j*(a*self.x+b*self.y)))
        return complex(value)

    def clock_c2(self, time, clock):
        return complex(np.sum(self.weights*clock.kernel(np.maximum(0,self.energy))*np.exp(-1j*self.energy*time)))

    def record(self):
        return {"family":self.smearing.family,"sigma_modular":self.smearing.sigma,
                "physical_smearing_rms_over_L":2*np.pi*self.smearing.sigma,
                "half_duration_modular":None if self.smearing.family=="gaussian" else self.smearing.half_duration,
                "frequency_radius":self.radius,"quadrature_nodes":len(self.energy),
                "mass_truncated":self.mass_truncated,"physical_variance_L1":self.physical_variance,
                "discarded_fraction_bound":self.discarded_fraction_bound}

@dataclass(frozen=True)
class ExponentialClock:
    mean: float
    def kernel(self,q):
        return np.exp(-np.maximum(0,np.asarray(q))/self.mean)
    def density(self,p):
        values=np.asarray(p)
        return np.where(values>=0,np.exp(-values/self.mean)/self.mean,0.)

def matched_clocks(mean, time_rms):
    W = 2*np.pi/(math.sqrt(3)*time_rms)
    w = 1/(math.sqrt(12)*time_rms)
    if mean < max(W/2,5*w):
        raise ValueError("Budget too small for the two matched shifted families")
    return {"compact":CompactClock(mean-W/2,W),
            "gamma5":SmoothClock(mean-5*w,w),
            "exponential_energy_only":ExponentialClock(mean)}

