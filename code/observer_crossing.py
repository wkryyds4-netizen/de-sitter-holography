"""Finite-mass dS2 observer crossed four-point kernel, arXiv:2411.08090v1.

A.14 and B.21--B.24, with equal real principal-series matter masses.
The returned values are complex; reality is a validation, never a projection.
No high-dimensional integral, coincident-label limit, or error certification
is supplied by this module.  HH values use the unnormalized f=sech(pi*s).
"""
from __future__ import annotations

from dataclasses import dataclass
from time import perf_counter
import mpmath


@dataclass
class CrossingResult:
    labels: tuple
    matter_nu: object
    raw_s: object
    gamma_product: object
    symmetric_s: object
    hh_cross: object
    g_terms: tuple
    cancellation_ratio: object
    requested_dps: int
    working_dps: int
    hyper_calls: int
    hyper_excess_errors: tuple
    elapsed_seconds: float


class CrossingEvaluator:
    """A private mpmath context avoids changing global precision.

    Inputs should be decimal strings for reproducible precision refinement.
    Only nondegenerate positive labels are supported.  The B.23 factors have
    removable poles on some coincident-label surfaces; we explicitly reject
    a=c or b=d instead of silently inserting zero or perturbing the input.
    """

    def __init__(self, dps=25, guard_digits=10, sum_method='r+s+e'):
        if dps < 15 or guard_digits < 0:
            raise ValueError('dps >= 15 and nonnegative guard_digits required')
        self.dps = int(dps)
        self.guard_digits = int(guard_digits)
        self.ctx = mpmath.mp.clone()
        self.ctx.dps = self.dps + self.guard_digits
        self.sum_method = sum_method
        self.hyper_calls = 0
        self.excess_errors = []

    def beta_like(self, u, v):
        """B.23; this symmetric Gamma ratio is NOT Euler beta."""
        c = self.ctx
        return c.sqrt(c.pi) * c.gammaprod(
            [(1-u)/2, (1-v)/2, (u+v-1)/2],
            [u/2, v/2, (2-u-v)/2])

    def hyper_f(self, h1, h2, h3, h4, h, k):
        """B.24: the 4F3 parameter excess is identically one."""
        c = self.ctx
        upper = [1+k-h1-h4, k+h1-h4, 2-h3-h4-h, 1+h3-h4-h]
        lower = [2-2*h4, 2+k-h2-h4-h, 1+k+h2-h4-h]
        self.excess_errors.append(abs(c.fsum(lower)-c.fsum(upper)-1))
        self.hyper_calls += 1
        return c.hyper(upper, lower, 1, strict=True, sum_method=self.sum_method)

    def g_terms(self, h1, h2, h3, h4, h, k):
        """The two terms in B.22, retained for cancellation diagnostics."""
        b = self.beta_like
        t1 = (b(2-h-h3-h4, 2*h4)
              * b(1+h-h1-h2, 1+h1-h4-k)
              * b(h2+h3-k, 1+h-h3-h4)
              * self.hyper_f(h1,h2,h3,1-h4,h,k))
        t2 = (b(h1+h4+k-1, 1-h1+h4-k)
              * b(1+h-h1-h2, h1+h4-k)
              * b(h2+h3-k, h-h3+h4)
              * self.hyper_f(h1,h2,h3,h4,h,k))
        return t1, t2

    def s_terms(self, r, nu):
        """B.21 with s_A=s_B=nu. Returns all four constituent terms."""
        c = self.ctx
        r1,r2,r3,r4 = r
        half = c.mpf('0.5')
        h1 = half + c.j*half*(-r2+r4-2*nu)
        h2 = half + c.j*half*(r1+r3)
        h3 = half + c.j*half*(-r2+r4+2*nu)
        h4 = half + c.j*half*(-r1-r3)
        h = half + c.j*half*(r1+r2-r3+r4)
        k = half + c.j*half*(-r1+r2+r3+r4)
        return (self.g_terms(h1,h2,h3,h4,h,k)
                + self.g_terms(h3,h4,h1,h2,1-h,1-k))

    def gamma_product(self, r, nu):
        """A.14's sixteen complex Gamma factors, via loggamma."""
        c = self.ctx
        quarter = c.mpf('0.25')
        half = c.mpf('0.5')
        log_terms = []
        for j, eta in enumerate((-1,-1,1,1)):
            x,y = r[j],r[(j+1)%4]
            for eps in (-1,1):
                log_terms.append(c.loggamma(quarter+c.j*half*(eta*nu+eps*(x+y))))
                log_terms.append(c.loggamma(quarter+c.j*half*(eps*nu+x-y)))
        return c.exp(c.fsum(log_terms))

    def gamma_product_expanded(self, r, nu):
        """Independent A.14 transcription, direct products (validation only)."""
        c = self.ctx
        a,b,z,d = r
        arguments = []
        for pair, matter_sign in [((a,b),-1),((b,z),-1),((z,d),1),((d,a),1)]:
            x,y = pair
            arguments.extend([
                (c.mpf('0.5')+c.j*(matter_sign*nu+x+y))/2,
                (c.mpf('0.5')+c.j*(matter_sign*nu-x-y))/2,
                (c.mpf('0.5')+c.j*(nu+x-y))/2,
                (c.mpf('0.5')+c.j*(-nu+x-y))/2,
            ])
        return c.fprod(c.gamma(z) for z in arguments)

    def evaluate(self, nu, a, b, z, d):
        c = self.ctx
        r = tuple(c.mpf(x) for x in (a,b,z,d))
        nu = c.mpf(nu)
        if nu <= 0 or any(x <= 0 for x in r):
            raise ValueError('strictly positive real principal-series labels required')
        if r[0] == r[2] or r[1] == r[3]:
            raise ValueError('coincident opposite labels need a removable-pole limit; unsupported')
        start = perf_counter()
        self.hyper_calls = 0
        self.excess_errors = []
        terms = self.s_terms(r,nu)
        term_sum = c.fsum(terms)
        raw_s = term_sum/c.pi**3
        gamma_product = self.gamma_product(r,nu)
        symmetric_s = raw_s*gamma_product/(4096*c.pi**11)
        a,b,z,d = r
        prefactor = a*c.sinh(c.pi*a)*b*c.tanh(c.pi*b)*z*c.sinh(c.pi*z)*d*c.tanh(c.pi*d)
        hh_cross = prefactor*symmetric_s
        cancellation = c.fsum(abs(t) for t in terms)/abs(term_sum) if term_sum else c.inf
        if not all(c.isfinite(x) for x in (raw_s,gamma_product,symmetric_s,hh_cross)):
            raise ArithmeticError('nonfinite crossing result; no real-part fallback is allowed')
        return CrossingResult(r,nu,raw_s,gamma_product,symmetric_s,hh_cross,
                              terms,cancellation,self.dps,c.dps,self.hyper_calls,
                              tuple(self.excess_errors),perf_counter()-start)

    def rho_hh(self, s):
        c = self.ctx
        s = c.mpf(s)
        return s*c.tanh(c.pi*s)/(4*c.pi**2*c.cosh(c.pi*s)**2)

    def conditional_cross(self, nu, a, b, z, d):
        """X_H(a,b,c,d)/rho_H(d), to multiply by normalized P_f(d)."""
        result = self.evaluate(nu,a,b,z,d)
        return result.hh_cross/self.rho_hh(d)
