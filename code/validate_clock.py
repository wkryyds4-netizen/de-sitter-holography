"""Independent quadrature, constrained primal and Fourier validation.

Paths are relative to this script. Finite grids are numerical evidence;
the continuum global-optimality proof remains necessary.
"""
from pathlib import Path
import json
import math
import platform
import warnings
import hashlib

import numpy as np
import scipy
from scipy.integrate import quad
from scipy.optimize import minimize, Bounds, LinearConstraint
from scipy.special import xlogy, logsumexp
from clock_model import Clock, CompactClock, SmoothClock, optimal_clock

ROOT = Path(__file__).resolve().parent
if ROOT.name == "code":
    ROOT = ROOT.parent
RECORDS = []


def record(name, passed, **details):
    RECORDS.append({"name": name, "passed": bool(passed), **details})


def close(name, got, expected, atol=2e-9, rtol=2e-8):
    record(name, math.isfinite(float(got)) and np.isclose(got, expected, atol=atol, rtol=rtol),
           got=float(got), expected=float(expected), atol=atol, rtol=rtol)


def integrate(fun, splits):
    return sum(quad(fun, a, b, epsabs=1e-11, epsrel=1e-10, limit=250)[0]
               for a, b in zip(splits[:-1], splits[1:]))


def density_integrals(obj, name, support):
    density = lambda p: float(obj.density(p))
    close(name + "/normalization", integrate(density, support), 1)
    close(name + "/mean", integrate(lambda p: p * density(p), support), obj.mean)
    def entropy(p):
        value = density(p)
        return value * (math.log(value) + p) if value > 0 else 0.0
    close(name + "/kl", integrate(entropy, support), obj.kl)
    Q = obj.threshold
    probes = [-1, 0, 0.2 * Q, Q, Q + 0.3, Q + 3]
    if isinstance(obj, CompactClock):
        probes += [obj.maximum_energy, obj.maximum_energy + 0.1]
    for i, q in enumerate(probes):
        left = max(0, q)
        endpoints = [left] + [z for z in support if z > left]
        value = integrate(density, endpoints) if len(endpoints) > 1 else 0.0
        close(name + f"/kernel-{i}", value, float(obj.kernel(q)))


def check_primal(Q, k, B, nlow):
    """Optimize bin probabilities, with reference-conditional density in each bin.

    Exact reference bin masses and conditional means make each feasible grid
    distribution continuum-feasible. Initialization does not use the model.
    """
    edges = np.r_[np.linspace(0, Q, nlow + 1),
                  np.linspace(Q, 16., 4 * nlow + 1)[1:]]
    lo, hi = edges[:-1], edges[1:]
    ref = np.exp(-lo) - np.exp(-hi)
    means = ((lo + 1) * np.exp(-lo) - (hi + 1) * np.exp(-hi)) / ref
    tail = (lo >= Q - 1e-12).astype(float)
    size, epsilon = len(ref), 1e-13
    # Generic discrete convex dual. This uses bin moments only, not the
    # continuum solution or its four analytic branches.
    def dual(z):
        eta, theta = z
        logits = np.log(ref)-eta*means+theta*tail
        logz = logsumexp(logits)
        probs = np.exp(logits-logz)
        return float(logz+eta*B-theta*k), np.array([B-means@probs, tail@probs-k])
    dual_result = minimize(dual, np.array([0.,0.]), jac=True, method="L-BFGS-B",
                           bounds=[(0,None),(0,None)],
                           options={"ftol":1e-15,"gtol":1e-12,"maxiter":1000,"maxls":50})
    eta, theta = dual_result.x
    logits = np.log(ref)-eta*means+theta*tail
    initial = np.exp(logits-logsumexp(logits))
    initial = np.maximum(initial,epsilon)
    initial /= np.sum(initial)
    objective = lambda x: float(np.sum(xlogy(x, x / ref)))
    jac = lambda x: np.log(np.maximum(x, 1e-300) / ref) + 1
    constraints = [LinearConstraint(np.ones((1, size)), 1, 1),
                   LinearConstraint(tail[None, :], k, np.inf),
                   LinearConstraint(means[None, :], -np.inf, B)]
    with warnings.catch_warnings(record=True) as observed:
        warnings.simplefilter("always")
        result = minimize(objective, initial, jac=jac, method="SLSQP",
                          bounds=Bounds(np.full(size, epsilon), np.ones(size)),
                          constraints=constraints,
                          options={"ftol": 2e-10, "maxiter": 1600})
    violation = max(abs(np.sum(result.x)-1), max(0, k-tail@result.x),
                    max(0, means@result.x-B))
    analytic = optimal_clock(Q, k, B)
    return {"bins": size, "success": bool(result.success), "message": str(result.message),
            "iterations": int(result.nit), "objective": float(result.fun),
            "analytic": analytic.kl, "gap": float(result.fun-analytic.kl),
            "mean": float(means@result.x), "tail": float(tail@result.x),
            "constraint_violation": float(violation),
            "dual_success":bool(dual_result.success), "dual_lower_bound":float(-dual_result.fun),
            "primal_dual_gap":float(result.fun+dual_result.fun),
            "warnings": sorted(set(str(w.message) for w in observed))}


def basic_and_kkt():
    cases = [(2,.1,2,"trace"), (2,.8,4,"tail_only"), (2,.02,.6,"mean_only"),
             (2,.8,2,"both"), (2,.1,.5,"both")]
    for i, (Q,k,B,regime) in enumerate(cases):
        obj = optimal_clock(Q,k,B)
        name = f"branches/{i}-{regime}"
        record(name+"/classification", obj.regime==regime, got=obj.regime, expected=regime)
        density_integrals(obj, name, [0,Q,math.inf])
        eta = obj.rate-1
        theta = math.log(obj.tail/(1-obj.tail))+math.log(math.expm1(obj.rate*Q))
        record(name+"/dual-feasibility", eta>=-1e-10 and theta>=-1e-10, eta=eta,theta=theta)
        close(name+"/mean-complementarity",eta*(B-obj.mean),0,atol=3e-9)
        close(name+"/tail-complementarity",theta*(obj.tail-k),0,atol=3e-9)
        p = np.array([.1,.4,Q-.1,Q+.2,Q+1,Q+3])
        stationarity = np.log(obj.density(p))+p+eta*p-theta*(p>=Q)
        close(name+"/stationarity",np.ptp(stationarity),0,atol=3e-9)
    tiny = Clock(1e-5,.73,.04,"manual")
    density_integrals(tiny,"tiny-threshold",[0,tiny.threshold,math.inf])
    for Q,k in [(2,.8),(.4,.2)]:
        for B in [0,k*Q*.99,k*Q]:
            try:
                optimal_clock(Q,k,B)
            except ValueError:
                record("infeasible",True,Q=Q,k=k,B=B)
            else:
                record("infeasible",False,Q=Q,k=k,B=B)
    Q,k = 2.,.8
    low_norm = quad(lambda p:math.exp(-p),0,Q)[0]
    b0 = (1-k)*quad(lambda p:p*math.exp(-p),0,Q)[0]/low_norm+k*(Q+1)
    for B in [b0,.6]:
        target_k = k if B==b0 else math.exp(-Q/B)
        left = optimal_clock(Q,target_k,B*(1-1e-7))
        right = optimal_clock(Q,target_k,B*(1+1e-7))
        close("boundary-continuity",left.kl,right.kl,atol=1e-6,rtol=1e-6)
    for B in [1-1e-8,1,1+1e-8]:
        close("trace-mean-boundary",optimal_clock(2,.01,B).kl,0,atol=1e-12)


def smooth_and_compact():
    for w in [.3,1.,3.]:
        obj = SmoothClock(2,w)
        name = f"gamma5/w={w}"
        density_integrals(obj,name,[0,2,math.inf])
        amp = lambda u:u*u*math.exp(-u/(2*w))/math.sqrt(24*w**5)
        der = lambda u:(2*u-u*u/(2*w))*math.exp(-u/(2*w))/math.sqrt(24*w**5)
        close(name+"/derivative-variance",quad(lambda u:der(u)**2,0,math.inf)[0],obj.time_variance)
        for tw in [0,.2,.5,1,3]:
            t = tw/w
            re = quad(lambda u:amp(u)*math.cos(t*u),0,math.inf,epsabs=1e-11,limit=400)[0]
            im = quad(lambda u:amp(u)*math.sin(t*u),0,math.inf,epsabs=1e-11,limit=400)[0]
            close(name+f"/fourier-{tw}",(re*re+im*im)/(2*math.pi),
                  float(obj.time_density(t)),atol=1e-10)
        close(name+"/time-normalization",quad(lambda t:float(obj.time_density(t)),-math.inf,math.inf)[0],1)
        close(name+"/time-variance",quad(lambda t:t*t*float(obj.time_density(t)),-math.inf,math.inf)[0],
              obj.time_variance)
    for W in [.7,2.,7.]:
        obj = CompactClock(2,W)
        name = f"compact/W={W}"
        density_integrals(obj,name,[0,2,2+W])
        amp = lambda u:math.sqrt(8/(3*W))*math.sin(math.pi*u/W)**2
        der = lambda u:math.sqrt(8/(3*W))*math.pi/W*math.sin(2*math.pi*u/W)
        close(name+"/derivative-variance",quad(lambda u:der(u)**2,0,W)[0],obj.time_variance)
        entropy = -quad(lambda u:xlogy(amp(u)**2,amp(u)**2),0,W)[0]
        close(name+"/entropy-closed-form",entropy,math.log(6*W)-7/3)
        for tw in [0,.4,2*math.pi,4*math.pi,7.5,35]:
            t = tw/W
            re = quad(lambda u:amp(u)*math.cos(t*u),0,W,epsabs=1e-11)[0]
            im = quad(lambda u:amp(u)*math.sin(t*u),0,W,epsabs=1e-11)[0]
            close(name+f"/fourier-{tw}",(re*re+im*im)/(2*math.pi),
                  float(obj.time_density(t)),atol=1e-11)
        # z=tW/(2pi). Each period is separately integrated. Envelope bounds tails.
        Z = 80
        norm = 2*sum(quad(lambda z:float(obj.time_density(2*math.pi*z/W))*2*math.pi/W,
                          j,j+1)[0] for j in range(Z))
        variance = 2*sum(quad(lambda z:(2*math.pi*z/W)**2*
                              float(obj.time_density(2*math.pi*z/W))*2*math.pi/W,
                              j,j+1)[0] for j in range(Z))
        tail_norm = 4/(15*math.pi**2*(1-Z**-2)**2*Z**5)
        tail_variance = 16/(9*W**2*(1-Z**-2)**2*Z**3)
        record(name+"/fourier-normalization-with-tail",abs(norm-1)<=tail_norm+1e-10,
               finite_integral=norm,tail_bound=tail_norm)
        record(name+"/fourier-variance-with-tail",abs(variance-obj.time_variance)<=tail_variance+1e-10,
               finite_integral=variance,expected=obj.time_variance,tail_bound=tail_variance)


def spectral_checks():
    clock = optimal_clock(3,.8)
    energies = np.array([2.8,2.9,-1.])
    proj_integral = quad(lambda p:float(clock.density(p))*float(np.all(p>=energies)),
                         0,max(energies),epsabs=1e-11)[0]
    proj_integral += quad(lambda p:float(clock.density(p)),max(energies),3)[0]
    proj_integral += quad(lambda p:float(clock.density(p)),3,math.inf)[0]
    close("n4/projector-max",proj_integral,float(clock.kernel(max(0,*energies))))
    sum_kernel = float(clock.kernel(np.sum(np.maximum(0,energies))))
    record("n4/max-is-not-sum",abs(proj_integral-sum_kernel)>.1,
           max_kernel=proj_integral,sum_kernel=sum_kernel)
    m = np.array([1.,-1.],dtype=complex)
    bulk,dressed = np.sum(m),np.sum(clock.kernel(np.array([0.,3.]))*m)
    record("relative-error-counterexample",abs(bulk)==0 and abs(dressed)>.1,
           bulk=float(bulk.real),dressed=float(dressed.real),L1=float(np.sum(abs(m))))
    rng = np.random.default_rng(261014175)
    E = rng.uniform(-2,3,(150,3))
    m = rng.normal(size=150)+1j*rng.normal(size=150)
    m /= np.sum(abs(m))
    K = clock.kernel(np.maximum(0,E.max(axis=1)))
    bound = (1-clock.tail)*float(np.sum(abs(m)))
    largest,l1_error = 0.,0.
    for scale in [0,1,10,1e2,1e4]:
        for _ in range(50):
            times = rng.normal(size=4)*scale
            evolved = m*np.exp(1j*(E@np.diff(times)))
            l1_error = max(l1_error,abs(np.sum(abs(evolved))-1))
            largest = max(largest,abs(np.sum((K-1)*evolved)))
    close("fixed-spectrum/L1",l1_error,0)
    record("fixed-spectrum/uniform-real-time-absolute-bound",largest<=bound+1e-12,
           largest_error=largest,bound=bound,samples=250,
           note="Numerical stress test; all-real-time uniformity follows from the proof.")
    E2 = np.vstack([E,[[0,7,0]]])
    m2 = np.r_[m*.9,.1j]
    q = np.maximum(0,E2.max(axis=1))
    tail_mass = float(np.sum(abs(m2[q>3])))
    bound2 = (1-clock.tail)*float(np.sum(abs(m2[q<=3])))+tail_mass
    err2 = abs(np.sum((clock.kernel(q)-1)*m2))
    record("spectral-tail-bound",err2<=bound2+1e-12,error=err2,bound=bound2,spectral_tail=tail_mass)


def main():
    basic_and_kkt()
    smooth_and_compact()
    spectral_checks()
    primal_results = {}
    for name,Q,k,B in [("trace",2,.1,2),("tail",2,.8,4),("mean",2,.02,.6),
                        ("both-high-k",2,.8,2),("both-low-k",2,.1,.5)]:
        rows = [check_primal(Q,k,B,n) for n in [6,12,24,48]]
        primal_results[name] = rows
        for row in rows:
            record("primal/"+name+f"/{row['bins']}",
                   row["dual_success"] and row["constraint_violation"]<1e-7 and row["gap"]>=-2e-7
                   and -2e-7<=row["primal_dual_gap"]<2e-7,
                   acceptance="Feasibility and primal-dual certificate; solver exit flag retained, not used as proof.",
                   **row)
        record("primal/"+name+"/refinement",
               rows[-1]["gap"]<=rows[0]["gap"]+2e-6 and rows[-1]["gap"]<.025,
               initial_gap=rows[0]["gap"],final_gap=rows[-1]["gap"],
               note="Truncated conditional-reference grid; continuum feasibility is exact.")
    failures = [x for x in RECORDS if not x["passed"]]
    output = {"status":"passed" if not failures else "failed",
              "code_sha256": {name: hashlib.sha256((Path(__file__).resolve().parent/name).read_bytes()).hexdigest()
                              for name in ("clock_model.py", "validate_clock.py")},
              "python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__,
              "checks":len(RECORDS),"failed":len(failures),"failures":failures,
              "primal":primal_results,"records":RECORDS,
              "coverage":["All four KL branches, including k<exp(-Q), boundaries and infeasibility",
                          "Independent quadrature: normalization, mean, KL and tail kernels",
                          "Independent SLSQP primal on reference-conditional grids with refinement and discrete dual certificate",
                          "KKT stationarity, multiplier signs and complementary slackness",
                          "Gamma5 and compact Fourier transforms, normalization and time variance",
                          "Compact entropy and hard energy support",
                          "n=4 max versus sum, cancellation, spectral tails and real-time phases"],
              "limitations":["Finite-grid numerical global checks supplement, not replace, the analytic proof.",
                             "Mathematical spectral probes are not reconstructed QFT observables.",
                             "No finite-G gravity, observer recoil, spatial localization or anti-scrambling amplitude."]}
    (ROOT/"results").mkdir(exist_ok=True)
    (ROOT/"results"/"validation.json").write_text(json.dumps(output,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({"status":output["status"],"checks":len(RECORDS),"failed":len(failures),
                      "failures":failures,"primal_final_gaps":{n:r[-1]["gap"] for n,r in primal_results.items()}},
                     ensure_ascii=False,indent=2))
    raise SystemExit(1 if failures else 0)


if __name__=="__main__":
    main()
