"""Regular real-label crossing kernel for arXiv:2411.08090v1.

Keeps observer_crossing.py unchanged.  Uses cyclic charts for a single
opposite-label coincidence and the analytic pole-subtracted derivative
formula at a=c,b=d.  Numerical derivatives/precision checks are not interval
certification.  All complex phases are retained.  HH prefactors always use
the ORIGINAL labels, even when the symmetric part is evaluated in a chart.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction
from time import perf_counter
from observer_crossing import CrossingEvaluator


@dataclass
class RegularCrossingResult:
    labels: tuple
    matter_nu: object
    chart_labels: tuple
    chart_rotation: int
    method: str
    symmetric_s: object
    hh_cross: object
    raw_s_chart: object
    gamma_product_chart: object
    cancellation_ratio: object
    residue_relative_error: object
    diagonal_a0: object
    diagonal_b0: object
    diagonal_difference_derivative: object
    requested_dps: int
    working_dps: int
    max_observed_formula_dps: int
    gap_guard_digits: int
    refinement_count: int
    hyper_calls: int
    max_hyper_excess_error: object
    reality_relative_error: object
    elapsed_seconds: float


class RegularCrossingEvaluator(CrossingEvaluator):
    """Pointwise API for nu>0 and strictly positive finite real a,b,c,d.

    evaluate(nu,a,b,c,d).hh_cross is the unnormalized HH spectral density.
    .symmetric_s is the D4-symmetric A.14 function, not the HH density.
    .conditional_cross divides the density by rho_H(d).

    Precision is increased for near-coincident opposite pairs and for
    detected cancellation. max_working_dps caps the outer formula context;
    mpmath.diff and special functions may use higher internal precision.
    max_observed_formula_dps records derivative callback precision.
    If the outer max_working_dps would be exceeded or the
    numerical reality/residue diagnostics fail, ArithmeticError is raised.
    No silent real projection or small-offset replacement is performed.
    """

    def __init__(self,dps=25,guard_digits=15,max_working_dps=180,sum_method='r+s+e'):
        super().__init__(dps,guard_digits,sum_method)
        self.max_working_dps = int(max_working_dps)
        self.max_observed_formula_dps = 0
        if self.max_working_dps < self.dps+self.guard_digits:
            raise ValueError('max_working_dps must cover requested precision and guards')

    def hyper_f(self,*args):
        self.max_observed_formula_dps=max(self.max_observed_formula_dps,self.ctx.dps)
        return super().hyper_f(*args)

    @staticmethod
    def _exact_scalar(value):
        """Exact input equality, before a finite precision parsing context."""
        if hasattr(value,'_mpf_'):
            sign,mantissa,exponent,_=value._mpf_
            return Fraction((-1 if sign else 1)*mantissa)*Fraction(2)**exponent
        return Fraction(value)

    def beta_stripped(self,u,v):
        """B.23 with only Gamma((u+v-1)/2) removed, analytically."""
        c = self.ctx
        return c.sqrt(c.pi)*c.gammaprod(
            [(1-u)/2,(1-v)/2], [u/2,v/2,(2-u-v)/2])

    def stripped_g_terms(self,h1,h2,h3,h4,h,k):
        """Both B.22 terms stripped of their common simple-pole Gamma."""
        b,bs = self.beta_like,self.beta_stripped
        first = (b(2-h-h3-h4,2*h4)
                 *bs(1+h-h1-h2,1+h1-h4-k)
                 *b(h2+h3-k,1+h-h3-h4)
                 *self.hyper_f(h1,h2,h3,1-h4,h,k))
        second = (b(h1+h4+k-1,1-h1+h4-k)
                  *b(1+h-h1-h2,h1+h4-k)
                  *bs(h2+h3-k,h-h3+h4)
                  *self.hyper_f(h1,h2,h3,h4,h,k))
        return first,second

    def stripped_branches(self,r,nu,return_terms=False):
        """pi^3 S = Gamma(i delta/2) A + Gamma(-i delta/2) B.

        delta=r1-r3. A and B remain analytic on the coincident surface.
        This method deliberately permits complex nearby labels for
        differentiation; the public evaluate API accepts only real ones.
        """
        c = self.ctx
        r1,r2,r3,r4 = r
        half = c.mpf('0.5')
        h1=half+c.j*half*(-r2+r4-2*nu)
        h2=half+c.j*half*(r1+r3)
        h3=half+c.j*half*(-r2+r4+2*nu)
        h4=half-c.j*half*(r1+r3)
        h=half+c.j*half*(r1+r2-r3+r4)
        k=half+c.j*half*(-r1+r2+r3+r4)
        first=self.stripped_g_terms(h1,h2,h3,h4,h,k)
        second=self.stripped_g_terms(h3,h4,h1,h2,1-h,1-k)
        if return_terms:
            return first+second
        return c.fsum(first),c.fsum(second)

    def _diagonal(self,r,nu):
        c=self.ctx
        x,y,_,_=r
        a0,b0=self.stripped_branches(r,nu)
        scale=max(abs(a0),abs(b0),c.mpf('1e-1000'))
        residue=abs(a0-b0)/scale
        # Centered path keeps r1+r3 fixed and r1-r3=delta.
        def difference(delta):
            a,b=self.stripped_branches((x+delta/2,y,x-delta/2,y),nu)
            return a-b
        derivative=c.diff(difference,c.zero,method='step',direction=0,addprec=15)
        derivative_part=2*derivative/c.j
        finite_part=-c.euler*(a0+b0)
        total=derivative_part+finite_part
        raw_s=total/c.pi**3
        cancellation=(abs(derivative_part)+abs(finite_part))/abs(total) if total else c.inf
        return raw_s,cancellation,residue,a0,b0,derivative

    def _ordinary(self,r,nu):
        c=self.ctx
        # Evaluate the singular Gamma from the original label difference.
        # In the unstripped B.21 form that difference is reconstructed from
        # several h parameters; roundoff there costs a SECOND set of gap
        # digits. Explicit extraction prevents that avoidable loss.
        delta=r[0]-r[2]
        stripped=self.stripped_branches(r,nu,return_terms=True)
        positive=c.gamma(c.j*delta/2)
        negative=c.gamma(-c.j*delta/2)
        terms=(positive*stripped[0],positive*stripped[1],
               negative*stripped[2],negative*stripped[3])
        total=c.fsum(terms)
        cancellation=c.fsum(abs(v) for v in terms)/abs(total) if total else c.inf
        return total/c.pi**3,cancellation,c.zero,None,None,None

    def evaluate(self,nu,a,b,z,d):
        c=self.ctx
        inputs=(nu,a,b,z,d)
        if any(isinstance(v,complex) for v in inputs):
            raise ValueError('public API requires real inputs')
        # Parse decimal inputs before selecting a lower working precision.
        # Otherwise very small nonzero gaps could be rounded to exact zero.
        with c.workdps(self.max_working_dps):
            original=tuple(c.mpf(v) for v in (a,b,z,d))
            nu=c.mpf(nu)
        if not all(c.isfinite(v) and v>0 for v in (nu,)+original):
            raise ValueError('strictly positive finite real labels and nu required')
        for i,j in ((0,2),(1,3)):
            if original[i]==original[j] and self._exact_scalar(inputs[i+1])!=self._exact_scalar(inputs[j+1]):
                raise ArithmeticError('distinct opposite inputs merged at max_working_dps; input gap exceeds precision budget')
        start=perf_counter()
        gap_ac=abs(original[0]-original[2])
        gap_bd=abs(original[1]-original[3])
        gap=max(gap_ac,gap_bd)
        diagonal=gap==0
        rotation=0 if diagonal or gap_ac>=gap_bd else 1
        chart=original[rotation:]+original[:rotation]
        scale=max(c.one,*original)
        gap_guard=0 if diagonal else max(0,int(c.ceil(-c.log10(gap/scale))))
        working=self.dps+self.guard_digits+gap_guard
        self.hyper_calls=0
        self.excess_errors=[]
        self.max_observed_formula_dps=0
        refinements=0
        while True:
            if working>self.max_working_dps:
                raise ArithmeticError(f'precision budget exceeded: need {working}, cap {self.max_working_dps}')
            with c.workdps(working):
                values=self._diagonal(chart,nu) if diagonal else self._ordinary(chart,nu)
                raw_s,cancellation,residue,a0,b0,derivative=values
                gamma_product=self.gamma_product(chart,nu)
                symmetric_s=raw_s*gamma_product/(4096*c.pi**11)
                # Critical: chart rotation does not act on the HH prefactor.
                aa,bb,cc,dd=original
                prefactor=(aa*c.sinh(c.pi*aa)*bb*c.tanh(c.pi*bb)
                           *cc*c.sinh(c.pi*cc)*dd*c.tanh(c.pi*dd))
                hh_cross=prefactor*symmetric_s
                if not all(c.isfinite(v) for v in (raw_s,symmetric_s,hh_cross,cancellation,residue)):
                    raise ArithmeticError('nonfinite regular crossing calculation')
                reality=abs(c.im(symmetric_s))/max(abs(symmetric_s),c.mpf('1e-1000'))
                loss=max(0,int(c.ceil(c.log10(max(c.one,cancellation)))))
                needed=self.dps+max(self.guard_digits,loss+8)
                target=c.power(10,-self.dps)
                if needed<=working and reality<=target and residue<=target:
                    result=RegularCrossingResult(
                        original,nu,chart,rotation,
                        'diagonal_derivative' if diagonal else ('cyclic_chart' if rotation else 'direct_chart'),
                        symmetric_s,hh_cross,raw_s,gamma_product,cancellation,residue,
                        a0,b0,derivative,self.dps,working,self.max_observed_formula_dps,gap_guard,refinements,
                        self.hyper_calls,max(self.excess_errors,default=c.zero),reality,
                        perf_counter()-start)
                    return result
                if refinements>=2:
                    raise ArithmeticError(f'regularity diagnostics failed: reality={reality}, residue={residue}')
                working=max(working+15,needed)
                refinements+=1

    def conditional_cross(self,nu,a,b,z,d):
        result=self.evaluate(nu,a,b,z,d)
        with self.ctx.workdps(result.working_dps):
            return result.hh_cross/self.rho_hh(d)
