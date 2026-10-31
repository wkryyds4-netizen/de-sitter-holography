"""High-precision fallback for an ill-conditioned individual Racah value.

This does not make the vectorized quadrature interval-certified. Keep
old raw grids unchanged; use this for explicit selected-point repairs.
"""
import mpmath as mp
from ds3_racah_fast import physical
from ds3_crossing_dictionary import physical_scalar


def checked_scalar(b,d,external,condition_limit=1e10):
    value,condition=physical(b,d,external)
    if float(condition)<=condition_limit:
        return complex(value),dict(method='complex128',condition=float(condition))
    # Preserve binary64 input rather than silently rounding onto a singular surface.
    values=[]
    for precision in [55,75]:
        with mp.workdps(precision):
            values.append(physical_scalar(mp.mpf(b),mp.mpf(d),[mp.mpf(x) for x in external]))
    with mp.workdps(75):
        discrepancy=abs(values[1]-values[0])
        relative=discrepancy/max(abs(values[1]),mp.mpf('1e-70'))
        if relative>mp.mpf('1e-25'):
            raise ArithmeticError('Precision refinement failed; do not use this point')
        return complex(values[-1]),dict(method='mpmath55_75',condition=float(condition),
               refinement_absolute=str(discrepancy),refinement_relative=str(relative),
               reference=[str(values[-1].real),str(values[-1].imag)])
