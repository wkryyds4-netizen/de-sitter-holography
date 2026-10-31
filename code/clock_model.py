"""Finite-clock spectral-kernel optimization, not a finite-G gravity simulator.

p and Q are dimensionless modular energies; logs are natural.
The density P is relative to dp; f=P/exp(-p) is relative to the II_1 trace.
"""
from dataclasses import dataclass, asdict
import math
import numpy as np
from scipy.optimize import brentq
from scipy.special import gammaincc, gammaln, digamma


def lower_mean(rate, threshold):
    """Mean of an exponential conditioned on 0 <= p < Q."""
    x = rate * threshold
    if x < 1e-3:
        return threshold * (0.5 - x / 12 + x**3 / 720 - x**5 / 30240)
    return 1 / rate - (threshold / math.expm1(x) if x < 700 else 0.0)


def split_mean(rate, threshold, tail):
    return (1 - tail) * lower_mean(rate, threshold) + tail * (threshold + 1 / rate)


def split_kl(rate, threshold, tail):
    """D(P || exp(-p)) for a two-piece conditioned exponential."""
    mean = split_mean(rate, threshold, tail)
    log_a0 = math.log1p(-tail) + math.log(rate) - math.log(-math.expm1(-rate * threshold))
    log_a1 = math.log(tail) + math.log(rate) + rate * threshold
    return (1 - tail) * log_a0 + tail * log_a1 + (1 - rate) * mean


@dataclass(frozen=True)
class Clock:
    threshold: float
    tail: float
    rate: float
    regime: str

    @property
    def mean(self):
        return split_mean(self.rate, self.threshold, self.tail)

    @property
    def kl(self):
        return split_kl(self.rate, self.threshold, self.tail)

    def density(self, p):
        values = np.asarray(p, dtype=float)
        result = np.zeros_like(values)
        low = (values >= 0) & (values < self.threshold)
        high = values >= self.threshold
        result[low] = ((1 - self.tail) * self.rate * np.exp(-self.rate * values[low])
                       / (-math.expm1(-self.rate * self.threshold)))
        result[high] = self.tail * self.rate * np.exp(-self.rate * (values[high] - self.threshold))
        return result

    def kernel(self, q):
        """Tail probability at max(0,q), evaluated stably on both pieces."""
        values = np.maximum(0, np.asarray(q, dtype=float))
        result = np.empty_like(values)
        low = values < self.threshold
        high = ~low
        r = self.rate
        result[low] = self.tail + (1 - self.tail) * (
            np.exp(-r * values[low]) * (-np.expm1(-r * (self.threshold - values[low])))
            / (-math.expm1(-r * self.threshold)))
        result[high] = self.tail * np.exp(-r * (values[high] - self.threshold))
        return result

    def record(self):
        return {**asdict(self), "mean": self.mean, "kl_nats": self.kl,
                "worst_kernel_error": 1 - self.tail}


def optimal_clock(threshold, min_tail, mean_cap=math.inf):
    """KL minimizer subject to K(Q)>=k and <p><=B.

    Strictly 0<k<1, Q>0. Densities (no atoms); B<=kQ is infeasible.
    Handles inactive, mean-only, tail-only and jointly active constraints.
    """
    Q, k, B = float(threshold), float(min_tail), float(mean_cap)
    if not Q > 0 or not 0 < k < 1 or not B > k * Q:
        raise ValueError("Require Q>0, 0<k<1 and B>kQ for a density.")
    rate_mean = max(1.0, 1.0 / B)
    tail_mean = math.exp(-rate_mean * Q)
    if tail_mean >= k:
        return Clock(Q, tail_mean, rate_mean, "trace" if rate_mean == 1 else "mean_only")
    if B >= split_mean(1.0, Q, k):
        return Clock(Q, k, 1.0, "tail_only")
    upper = 2.0
    while split_mean(upper, Q, k) > B:
        upper *= 2
    rate = brentq(lambda z: split_mean(z, Q, k) - B, 1.0, upper,
                  xtol=1e-13, rtol=1e-13)
    return Clock(Q, k, rate, "both")


def exponential_baseline(threshold, min_tail):
    """Least KL in the original unshifted heavy exponential family 0<alpha<=1."""
    alpha = min(1.0, -math.log(min_tail) / threshold)
    return {"alpha": alpha, "mean": 1 / alpha,
            "kl_nats": math.log(alpha) + 1 / alpha - 1,
            "tail": math.exp(-alpha * threshold)}


def shifted_exponential(threshold):
    """k=1 limit: normal but non-faithful state, exact only on q<=Q."""
    return {"mean": threshold + 1, "kl_nats": threshold, "tail": 1.0}


@dataclass(frozen=True)
class SmoothClock:
    """Shifted Gamma(shape=5) pure real wave amplitude, with Fourier time POVM.

    Exact kernel recovery on q<=Q; finite timing variance. This does not
    impose locality, coupled gravity, a hard upper energy cutoff or faithfulness.
    """
    threshold: float
    scale: float

    def __post_init__(self):
        if self.threshold < 0 or self.scale <= 0:
            raise ValueError("Require Q>=0, w>0")

    @property
    def mean(self):
        return self.threshold + 5 * self.scale

    @property
    def time_variance(self):
        return 1 / (12 * self.scale**2)

    @property
    def kl(self):
        h = 5 + math.log(self.scale) + gammaln(5) - 4 * digamma(5)
        return float(self.mean - h)

    def density(self, p):
        u = np.asarray(p, dtype=float) - self.threshold
        result = np.zeros_like(u)
        mask = u > 0
        result[mask] = np.exp(4 * np.log(u[mask]) - u[mask] / self.scale
                              - gammaln(5) - 5 * math.log(self.scale))
        return result

    def kernel(self, q):
        return gammaincc(5, np.maximum(0, np.asarray(q, dtype=float) - self.threshold) / self.scale)

    def time_density(self, t):
        """|Fourier(sqrt(P))|^2 on the full-line dilation, integral = 1."""
        t = np.asarray(t, dtype=float)
        return 1 / (12 * np.pi * self.scale**5 * (1 / (4 * self.scale**2) + t*t)**3)

    def record(self):
        return {**asdict(self), "mean": self.mean, "kl_nats": self.kl,
                "time_rms": math.sqrt(self.time_variance),
                "worst_kernel_error": 0.0}


@dataclass(frozen=True)
class CompactClock:
    """Smooth compact-support amplitude sqrt(8/(3W))*sin^2(pi*(p-Q)/W).

    The full-line zero extension is H^2. The timing quantity is a covariant
    Fourier POVM variance, not a self-adjoint canonical time on p>=0.
    """
    threshold: float
    width: float

    def __post_init__(self):
        if self.threshold < 0 or self.width <= 0:
            raise ValueError("Require Q>=0, W>0")

    @property
    def mean(self):
        return self.threshold + self.width / 2

    @property
    def maximum_energy(self):
        return self.threshold + self.width

    @property
    def time_variance(self):
        return 4 * math.pi**2 / (3 * self.width**2)

    @property
    def kl(self):
        return self.mean - math.log(6 * self.width) + 7 / 3

    def density(self, p):
        x = (np.asarray(p, dtype=float) - self.threshold) / self.width
        return np.where((x > 0) & (x < 1), 8 / (3 * self.width) * np.sin(np.pi * x)**4, 0.0)

    def kernel(self, q):
        x = np.clip((np.asarray(q, dtype=float) - self.threshold) / self.width, 0, 1)
        cdf = x - 2 / (3 * np.pi) * np.sin(2 * np.pi * x) + 1 / (12 * np.pi) * np.sin(4 * np.pi * x)
        return np.clip(1 - cdf, 0, 1)

    def time_density(self, t):
        z = np.asarray(t, dtype=float) * self.width / (2 * np.pi)
        bracket = 0.5 * np.sinc(z) + 0.25 * (np.sinc(z - 1) + np.sinc(z + 1))
        return 4 * self.width / (3 * np.pi) * bracket**2

    def record(self):
        return {**asdict(self), "mean": self.mean, "kl_nats": self.kl,
                "max_energy": self.maximum_energy,
                "time_rms": math.sqrt(self.time_variance),
                "worst_kernel_error": 0.0}
