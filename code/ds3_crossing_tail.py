"""Conservative absolute bound on excluded Euclidean puncture regions."""
import mpmath as mp


def puncture_tail_bound(epsilon):
    """Raw-6j tail in all six chart regions r<epsilon.

    All dimensions have real part one. The shadow integral is bounded by
    I(z)=int d2x/(|x||x-z||x-1|). This is NOT a quadrature error bound.
    """
    ep = mp.mpf(epsilon)
    delta = ep/(1-ep)
    if not 0 < delta <= mp.mpf('.25'):
        raise ValueError('Require 0 < epsilon <= 1/5')
    v = 12+mp.log(1/delta)
    c = 128*mp.pi**3/3
    return c*(2*delta*(v*v+2*v+2)+delta**2/2*(v*v+v+mp.mpf('.5')))
