"""Pure endpoint-module helpers, epsilon = -log(q) = 4 lambda.

No scientific output is written on import. B/P/E use Lin--Stanford (48),
with sectorwise rescaling (147); Btr = (K_L-K_R)/2.
"""
from functools import lru_cache
import math
import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import expm_multiply
from scipy.special import gammaln
from scipy.stats import poisson


def qinteger(epsilon, n):
    return -np.expm1(-epsilon * np.asarray(n)) / (-np.expm1(-epsilon))


def heat_vector(epsilon, m, cutoff, method="expm"):
    """Return alternating normalized exp(-sqrt(m) T_q)|0>.

    Store the scaled positive vector exp(sqrt(m) T_q - m)|0>.
    Taylor uses positivity and an explicit finite-matrix operator-norm
    remainder bound, independent of expm_multiply's implementation.
    """
    n = np.arange(cutoff + 1)
    off = np.sqrt(qinteger(epsilon, n[1:]))
    T = diags([off, off], [-1, 1], shape=(cutoff+1, cutoff+1), format="csr")
    t = math.sqrt(m)
    start = np.zeros(cutoff+1)
    start[0] = 1
    meta = {}
    if method == "expm":
        # Keep the matrix entrywise positive during the action. Putting -m I
        # inside expm_multiply caused premature small-tail truncation in this
        # test (documented diagnostic). m <= 512 keeps exp(m) representable.
        x = expm_multiply(t*T, start, traceA=0)*math.exp(-m)
    elif method == "taylor":
        radius = 2*t*math.sqrt(cutoff)
        order = max(1, int(radius))
        while radius-m+poisson.logsf(order, radius) > -85:
            order += 1
        term = start*math.exp(-m)
        x = term.copy()
        for k in range(1, order+1):
            term = (t/k)*(T @ term)
            x += term
        meta = dict(order=order, scaled_vector_tail_bound=math.exp(
            radius-m+poisson.logsf(order, radius)))
    else:
        raise ValueError(method)
    norm = np.linalg.norm(x)
    return ((-1.)**n)*x/norm, dict(
        scaled_norm_squared=float(norm**2), minimum_positive_coefficient=float(x.min()), **meta)


def coherent(m, cutoff):
    n = np.arange(cutoff+1)
    return (-1.)**n * np.exp((-m+n*np.log(m)-gammaln(n+1))/2)


def endpoint_actions(epsilon, d, psi):
    """e_n and orthogonal f_n coordinates, including Btr's last raising row."""
    psi = np.r_[psi, 0.]
    n = np.arange(len(psi))
    qn = qinteger(epsilon, n)
    r = np.exp(-epsilon*d)
    omr = -np.expm1(-epsilon*d)
    omr2 = -np.expm1(-2*epsilon*d)
    c = np.exp(-epsilon*(n+d)/2)
    root = np.sqrt(-np.expm1(-epsilon*(n+d)))
    g = np.sqrt(qn*omr2)
    B = np.array([qn*omr2/(2*c*root)*psi, -r*g/(2*c*root)*psi])
    P = np.array([np.zeros_like(psi, dtype=complex), -1j*g/(2*root)*psi])
    E = np.array([(qinteger(epsilon,d)-omr**2*qn/2)/c*psi, -r*g/(2*c)*psi])
    prev = np.r_[0., psi[:-1]]
    nxt = np.r_[psi[1:], 0.]
    Q = np.array([omr/(2*np.sqrt(epsilon))*(np.sqrt(qn)*prev
                      +np.sqrt(qinteger(epsilon,n+1))*nxt),
                  -np.sqrt(omr2)/(2*np.sqrt(epsilon))*prev])
    state = np.array([psi, np.zeros_like(psi)])
    return dict(B=B, P=P, E=E, Btr=Q, state=state)


def observables(epsilon, d, psi):
    a = endpoint_actions(epsilon, d, psi)
    out = {}
    for name in ("B", "P", "E", "Btr"):
        mean = np.vdot(a["state"], a[name])
        second = np.vdot(a[name], a[name]).real
        out[name] = dict(mean_real=float(mean.real), mean_imag=float(mean.imag),
                         second=float(second), variance=float(second-abs(mean)**2))
    residual = np.linalg.norm(a["B"]+a["Btr"])
    out.update(boost_sum_norm=float(residual),
               boost_difference_norm=float(np.linalg.norm(a["B"]-a["Btr"])),
               boost_sum_cancellation=float((np.linalg.norm(a["B"])
                       +np.linalg.norm(a["Btr"]))/residual),
               E_eigenvector_defect=float(np.linalg.norm(a["E"]-d*a["state"])),
               lowering_defect=float(np.linalg.norm(a["B"]+1j*a["P"])),
               BP_commutator_expectation_imag=float(2*np.vdot(a["B"],a["P"]).imag))
    return out


def gram_factory(q, r):
    """Exact rational chord Gram via annihilation of a bra creation."""
    def qi(n): return (1-q**n)/(1-q)
    @lru_cache(None)
    def gram(u, v):
        if sum(u) != sum(v): return q*0
        if sum(u) == 0: return q*0+1
        l, k = v
        if u[0]:
            a=(u[0]-1,u[1])
            return ((qi(l)*gram(a,(l-1,k)) if l else 0)
                    +(r*q**l*qi(k)*gram(a,(l,k-1)) if k else 0))
        a=(0,u[1]-1)
        return ((qi(k)*gram(a,(l,k-1)) if k else 0)
                +(r*q**k*qi(l)*gram(a,(l-1,k)) if l else 0))
    return gram
