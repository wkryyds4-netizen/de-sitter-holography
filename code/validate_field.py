"""Stage-two independent field validation via adaptive convolution/CDF quadrature."""
from pathlib import Path
import json, math, cmath, hashlib, platform, time
from functools import lru_cache
import numpy as np
import scipy
from scipy.integrate import quad, quad_vec
from field_model import Smearing, Spectrum, physical_rho, matched_clocks

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/"results/stage2"
RECORDS=[]
WARNINGS=[]

def record(name, passed, **kwargs):
    RECORDS.append({"name":name,"passed":bool(passed),**kwargs})

def close(name,got,expected,atol=2e-8,rtol=2e-8,**kwargs):
    error=abs(got-expected)
    record(name,error<=atol+rtol*abs(expected),got=float(got),expected=float(expected),
           absolute_error=float(error),atol=atol,rtol=rtol,**kwargs)

def cclose(name,got,expected,atol=2e-8,rtol=2e-8):
    record(name,abs(got-expected)<=atol+rtol*abs(expected),
           got=[got.real,got.imag],expected=[expected.real,expected.imag],
           absolute_error=abs(got-expected),atol=atol,rtol=rtol)

def sinc(x):
    return math.sin(x)/x if abs(x)>1e-12 else 1-x*x/6

class ReferenceSpectrum:
    """Independent scalar density and adaptive, piecewise-integrated CDF."""
    def __init__(self,family,sigma,radius):
        self.family,self.sigma,self.R=family,sigma,radius
        self.T=sigma/math.sqrt(1/3-2/math.pi**2)
        self.step=math.pi/self.T if family=="cosine" else min(2,1/sigma)
        self.edges=np.unique(np.r_[0,np.arange(self.step,radius,self.step),radius])
        segments=[quad(self.positive,a,b,epsabs=1e-11,epsrel=2e-11)[0]
                  for a,b in zip(self.edges[:-1],self.edges[1:])]
        self.prefix=np.r_[0,np.cumsum(segments)]
        self.negative=quad(lambda x:math.exp(-x)*self.positive(x),0,radius,
                           epsabs=1e-11,epsrel=2e-11,limit=300)[0]
        self.Z=float(self.prefix[-1]+self.negative)
        self.F0=self.negative/self.Z
    def hat(self,E):
        if self.family=="gaussian":
            return math.exp(-.5*(self.sigma*E)**2)
        z=E*self.T
        if abs(abs(z)-math.pi)<.05 or abs(z)<.05:
            return sinc(z)+.5*(sinc(z-math.pi)+sinc(z+math.pi))
        return math.pi**2*math.sin(z)/(z*(math.pi**2-z*z))
    def positive(self,E):
        if abs(E)<1e-6:
            shape=1+E/2+E*E/12
        else:
            shape=E/(-math.expm1(-E))
        return shape*self.hat(E)**2
    def mu(self,E):
        if abs(E)>self.R:
            return 0.
        return self.positive(abs(E))*(math.exp(E) if E<0 else 1)/self.Z
    @lru_cache(maxsize=30000)
    def F(self,q):
        if q<=-self.R:return 0.
        if q>=self.R:return 1.
        if q<0:
            return quad(lambda u:math.exp(-u)*self.positive(u),-q,self.R,
                        epsabs=2e-11,epsrel=2e-10,limit=300)[0]/self.Z
        index=min(np.searchsorted(self.edges,q,side="right")-1,len(self.edges)-2)
        local=quad(self.positive,self.edges[index],q,epsabs=2e-11,epsrel=2e-10)[0]
        return (self.negative+self.prefix[index]+local)/self.Z
    @lru_cache(maxsize=30000)
    def q_densities(self,q):
        """q>0 densities of max(0,X,Y), X_++Y_+, (X+Y_+)_+."""
        u=self.mu(q)
        lo,hi=max(0,q-self.R),min(q,self.R)
        if lo<hi:
            mid=(lo+hi)/2
            H=(quad(lambda x:self.mu(x)*self.mu(q-x),lo,mid,
                    epsabs=2e-12,epsrel=2e-9,limit=300)[0]
               +quad(lambda x:self.mu(x)*self.mu(q-x),mid,hi,
                     epsabs=2e-12,epsrel=2e-9,limit=300)[0])
        else:H=0.
        if q<self.R:
            limit=self.R-q
            cut=min(12.,limit)
            Rterm=quad(lambda v:self.mu(q+v)*self.mu(-v),0,cut,
                       epsabs=2e-12,epsrel=2e-9,limit=200)[0]
            if limit>cut:
                Rterm+=quad(lambda v:self.mu(q+v)*self.mu(-v),cut,limit,
                            epsabs=2e-12,epsrel=2e-9,limit=200)[0]
        else:Rterm=0.
        return np.array([2*self.F(q)*u,2*self.F0*u+H,self.F0*u+H+Rterm])
    def integrate_clocks(self,clocks):
        points=np.unique(np.r_[self.edges,self.edges+self.R,
                               np.arange(0,min(32,2*self.R),.5),
                               [c.threshold for c in clocks if hasattr(c,"threshold")],
                               [c.maximum_energy for c in clocks if hasattr(c,"maximum_energy")]])
        points=points[(points>0)&(points<2*self.R)]
        def integrand(q):
            deficits=np.array([1-float(c.kernel(q)) for c in clocks])
            return np.r_[self.mu(q),self.q_densities(q)][:,None]*deficits[None,:]
        result,error,info=quad_vec(integrand,0,2*self.R,epsabs=2e-9,epsrel=2e-9,
                                   points=points,limit=10000,full_output=True)
        return result,float(error),{"evaluations":int(info.neval),"success":bool(info.success),
                                    "message":info.message}
    def Js(self,p):
        Fp=self.F(p)
        J1=Fp*Fp
        J2=self.F0*Fp+quad(lambda x:self.mu(x)*self.F(p-x),0,min(p,self.R),
                         epsabs=2e-10,epsrel=2e-9,limit=400)[0]
        upper=min(p,self.R)
        J3=quad(lambda x:self.mu(x)*self.F(p-x),-self.R,0,
                epsabs=2e-10,epsrel=2e-9,limit=400)[0]
        if upper>0:
            J3+=quad(lambda x:self.mu(x)*self.F(p-x),0,upper,
                     epsabs=2e-10,epsrel=2e-9,limit=400)[0]
        return np.array([J1,J2,J3])

def independent_field_basics(ref,model,meta):
    name=f"{ref.family}/{ref.sigma}"
    close(name+"/normalization",meta["mass_truncated"],ref.Z,atol=2e-9,rtol=2e-10)
    close(name+"/physical-variance",meta["physical_variance_L1"],ref.Z/(16*math.pi**4),
          atol=2e-11,rtol=2e-10)
    support=math.inf if ref.family=="gaussian" else ref.T
    close(name+"/time-smear-mass",2*quad(lambda t:float(model.time_function(t)),0,support)[0],1)
    close(name+"/time-smear-rms2",2*quad(lambda t:t*t*float(model.time_function(t)),0,support)[0],
          ref.sigma**2)
    for E in [0,.7,3,math.pi/ref.T,2*math.pi/ref.T]:
        direct=2*quad(lambda s:float(model.time_function(s))*math.cos(E*s),0,support,
                     epsabs=1e-11,limit=400)[0]
        close(name+f"/direct-smear-Fourier/{E}",direct,float(model.fourier(E)),atol=3e-10)
    for E in [.1,1,5,20]:
        close(name+f"/KMS-spectrum/{E}",float(model.density(-E)),
              math.exp(-E)*float(model.density(E)),atol=1e-15)
    # Independent Euclidean-shifted Wightman convolution avoids coincidence poles.
    # W(s-i theta)=-(16 pi²)^-1 csch²(pi(s-i theta)) for radius L=1.
    theta=.5
    def wightman(s):
        return -1/(16*math.pi**2*cmath.sinh(math.pi*complex(s,-theta))**2)
    def auto(v):
        if ref.family=="gaussian":
            return math.exp(-v*v/(4*ref.sigma**2))/(2*math.sqrt(math.pi)*ref.sigma)
        low=max(-ref.T,-ref.T-v); high=min(ref.T,ref.T-v)
        if high<=low:return 0.
        return quad(lambda u:float(model.time_function(u))*float(model.time_function(u+v)),
                    low,high,epsabs=1e-11)[0]
    extent=12*ref.sigma if ref.family=="gaussian" else 2*ref.T
    s=.17
    temporal=complex(quad(lambda v:auto(v)*wightman(s+v).real,-extent,extent,
                          epsabs=1e-11,limit=300)[0],
                     quad(lambda v:auto(v)*wightman(s+v).imag,-extent,extent,
                          epsabs=1e-11,limit=300)[0])
    # At theta=1/2 the negative and positive frequency contributions combine.
    spectral=quad(lambda E:2*ref.positive(E)*math.exp(-E/2)*math.cos(E*s)/(16*math.pi**4),
                  0,ref.R,epsabs=1e-11,limit=500)[0]
    cclose(name+"/Wightman-time-vs-spectral",temporal,complex(spectral),atol=2e-10)
    for t in [.0,.4,1.2]:
        def integral(sign,thermal=False):
            # Original W(t), or W(t-i). Explicit two signed-energy pieces.
            def real(E):
                h=ref.positive(E)/ref.Z
                return h*(1+math.exp(-E))*math.cos(E*t)
            def imag(E):
                h=ref.positive(E)/ref.Z
                return sign*h*(1-math.exp(-E))*math.sin(E*t)
            return complex(quad(real,0,ref.R,epsabs=2e-10,limit=600)[0],
                           quad(imag,0,ref.R,epsabs=2e-10,limit=600)[0])
        Wminus=integral(-1)
        # Evaluate the model density directly on both energy signs and the
        # imaginary-time boundary, independently of the split reference formula.
        Wshift=complex(
            quad(lambda E:float(model.density(E))*math.exp(-E)*math.cos(E*t)/ref.Z,
                 -ref.R,ref.R,points=[0],epsabs=1e-10,limit=1200)[0],
            quad(lambda E:-float(model.density(E))*math.exp(-E)*math.sin(E*t)/ref.Z,
                 -ref.R,ref.R,points=[0],epsabs=1e-10,limit=1200)[0])
        cclose(name+f"/KMS-Fourier/{t}",Wshift,Wminus.conjugate(),atol=2e-9)

def finite_mode_check():
    omega=1.
    N=48
    weights=np.exp(-omega*np.arange(N)); weights/=weights.sum()
    vacuum=np.diag(np.sqrt(weights))
    annih=np.diag(np.sqrt(np.arange(1,N)),1)
    field=annih+annih.T
    Z=float(np.sum(weights*np.diag(field@field)))
    field=field/math.sqrt(Z)
    energy=omega*(np.arange(N)[:,None]-np.arange(N)[None,:])
    cutoff=1.5*omega
    projector=(energy<=cutoff)
    def operator(t):
        return field*np.exp(1j*energy*t)
    def product(times):
        state=vacuum.astype(complex)
        for j,t in enumerate(reversed(times)):
            state=operator(t)@state
            if j<len(times)-1:state=state*projector
        return complex(np.vdot(vacuum,state))
    a=1/(1+math.exp(-omega)); b=1-a
    close("finite-mode/KMS-ratio",b/a,math.exp(-omega),atol=1e-14)
    close("finite-mode/C2-perfect",product([0,0]).real,1,atol=1e-12)
    close("finite-mode/C4-nonGaussian",product([0,0,0,0]).real,3-2*a*a,atol=1e-12)
    for t in [0,.7,2.]:
        total=0j
        for x,px in [(omega,a),(-omega,b)]:
            for y,py in [(omega,a),(-omega,b)]:
                q=(max(0,x,y),max(x,0)+max(y,0),max(0,x+max(y,0)))
                phases=(np.exp(-1j*t*(x+y)),1,np.exp(1j*t*(y-x)))
                total+=px*py*sum(phase for qi,phase in zip(q,phases) if qi<=cutoff)
        cclose(f"finite-mode/OTO-matrix-vs-Wick-measure/{t}",product([t,0,t,0]),total,atol=2e-12)

def time_tail_checks(cases):
    details={}
    for family in ("compact","gamma5"):
        row=next(r for r in cases if r["clock_family"]==family)
        clock=matched_clocks(row["clock_mean"],.5)[family]
        saved=row["clock"]["time_tail_abs_gt_1"]
        true=1-2*quad(lambda t:float(clock.time_density(t)),0,1,epsabs=1e-12)[0]
        if family=="compact":
            z=200*clock.width/(2*math.pi)
            bound=4/(15*math.pi**2*(1-z**-2)**2*z**5)
        else:
            bound=1/(30*math.pi*clock.scale**5*200**5)
        close(f"clock-tail/{family}/finite-vs-complement",saved,true,
              atol=2e-9+bound,rtol=0)
        details[family]={"saved_integral_1_to_200":saved,"independent_total_tail":true,
                         "absolute_difference":abs(saved-true),"omitted_gt_200_bound":bound}
    return details

def main():
    start=time.perf_counter()
    study=json.loads((OUT/"study.json").read_text())
    cases=study["cases"]
    comparisons=[]; tails=[]; refs={}
    # Check the physical prefactor without field or clock models.
    for theta in [.2,.5,.8]:
        integral=quad(lambda E: (E/(-math.expm1(-E)) if E>1e-8 else 1)*
                      (math.exp(-theta*E)+math.exp(-(1-theta)*E))/(16*math.pi**4),
                      0,math.inf,epsabs=1e-12)[0]
        close(f"physical-spectrum/Euclidean-{theta}",integral,
              1/(16*math.pi**2*math.sin(math.pi*theta)**2),atol=1e-12)
    for E in [-3,0,1,5]:
        shape=(E/(-math.expm1(-E))) if E else 1
        close(f"physical-spectrum/prefactor-{E}",float(physical_rho(E,2.3)),
              shape/(16*math.pi**4*2.3**2),atol=1e-15)
    for meta in study["spectra"]:
        family,sigma=meta["family"],meta["sigma_modular"]
        model=Smearing(family,sigma)
        ref=ReferenceSpectrum(family,sigma,meta["frequency_radius"])
        refs[(family,sigma)]=ref
        independent_field_basics(ref,model,meta)
        selected=[r for r in cases if r["field_family"]==family and r["field_sigma"]==sigma]
        clocks=[matched_clocks(r["clock_mean"],.5)[r["clock_family"]] for r in selected]
        values,outer_error,info=ref.integrate_clocks(clocks)
        record(f"adaptive-reference/{family}/{sigma}",info["success"] and outer_error<1e-7,
               outer_error_estimate=outer_error,**info)
        for j,row in enumerate(selected):
            ref2=float(values[0,j]); pieces=values[1:,j]; ref4=float(pieces.mean())
            diff2=row["sup_C2_absolute"]-ref2
            diff4=row["sup_C4_absolute_over_3"]-ref4
            details={"field_family":family,"field_sigma":sigma,"clock_family":row["clock_family"],
                     "clock_mean":row["clock_mean"],"reference_C2":ref2,"reference_C4_over_3":ref4,
                     "saved_C2":row["sup_C2_absolute"],"saved_C4_over_3":row["sup_C4_absolute_over_3"],
                     "C2_signed_difference":diff2,"C4_signed_difference":diff4,
                     "reference_pairing_deficits":pieces.tolist(),
                     "saved_pairing_deficits":row["C4_pairing_deficits"],
                     "outer_error_estimate":outer_error,
                     "spectral_cutoff_note":"Independent adaptive calculation at the same R; omitted tails bounded in source study."}
            comparisons.append(details)
            # Explicit absolute tolerance; no relative-accuracy claim for tiny tails.
            record("72-case/C2",abs(diff2)<=2e-6,absolute_tolerance=2e-6,**details)
            record("72-case/C4",abs(diff4)<=2e-6,absolute_tolerance=2e-6,**details)
        if sigma==.3:
            clock=matched_clocks(10,.5)["compact"]
            Q=clock.threshold
            Js=ref.Js(Q); Tpieces=1-Js
            t2=1-ref.F(Q); t4=float(Tpieces.mean())
            chosen=next(c for c in comparisons if c["field_family"]==family and c["field_sigma"]==sigma
                        and c["clock_family"]=="compact" and c["clock_mean"]==10)
            # Compare the repaired production CDF tails with independent adaptive CDF integrals.
            gridtails=Spectrum(model,order=20,periods=50).tails(Q)
            entry={"field_family":family,"field_sigma":sigma,"clock":"compact","mean":10,"Q":Q,
                   "delta_Q":0.,"CDF_T2":t2,"CDF_T4_over_3":t4,"CDF_T4_pieces":Tpieces.tolist(),
                   "actual_C2_error":chosen["reference_C2"],"actual_C4_error_over_3":chosen["reference_C4_over_3"],
                   "window_C2_bound":t2,"window_C4_bound":t4,"production_CDF_tail":gridtails,
                   "production_CDF_C2_signed_error":gridtails["T2_normalized"]-t2,
                   "production_CDF_C4_signed_error":gridtails["T4_normalized"]-t4}
            tails.append(entry)
            close("production-CDF-T2/"+family,gridtails["T2_normalized"],t2,atol=2e-7)
            close("production-CDF-T4/"+family,gridtails["T4_normalized"],t4,atol=2e-7)
            for i in range(3):
                close(f"production-CDF-pair/{family}/{i}",gridtails["T4_pieces"][i],Tpieces[i],atol=2e-7)
            record("CDF-tail-bound/"+family,chosen["reference_C2"]<=t2+2e-8 and
                   chosen["reference_C4_over_3"]<=t4+2e-8,**entry)
            # Independent derivative-density integral matches direct 2D CDF.
            intpieces=np.array([quad(lambda q,i=i:ref.q_densities(q)[i],Q,2*ref.R,
                                    epsabs=2e-8,limit=500)[0] for i in range(3)])
            for i in range(3):
                close(f"CDF-vs-convolution-tail/{family}/{i}",intpieces[i],Tpieces[i],atol=2e-7)
        print(f"Independent adaptive validation completed {family} sigma={sigma}",flush=True)
    # Independent second route: average the direct CDFs over the clock
    # density, rather than integrating the q convolution densities.
    cdf_crosschecks=[]
    for family,sigma,mean,kind in [("cosine",.1,5.,"compact"),
                                   ("gaussian",.3,10.,"compact"),
                                   ("cosine",.3,10.,"compact"),
                                   ("cosine",.1,5.,"gamma5")]:
        ref=refs[(family,sigma)]
        clock=matched_clocks(mean,.5)[kind]
        low=clock.threshold
        high=clock.maximum_energy if kind=="compact" else min(2*ref.R,low+80*clock.scale)
        direct,error,info=quad_vec(lambda p:float(clock.density(p))*(1-ref.Js(p)),low,high,
                                  epsabs=3e-8,epsrel=3e-8,full_output=True,limit=1000)
        row=next(r for r in comparisons if r["field_family"]==family and r["field_sigma"]==sigma
                 and r["clock_mean"]==mean and r["clock_family"]==kind)
        residual=np.asarray(row["reference_pairing_deficits"])-direct
        entry={"family":family,"sigma":sigma,"mean":mean,"clock":kind,
               "CDF_clock_average":direct.tolist(),"convolution_deficits":row["reference_pairing_deficits"],
               "max_residual":float(max(abs(residual))),"quadrature_error_estimate":float(error)}
        cdf_crosschecks.append(entry)
        record("clock-CDF-vs-q-convolution",info.success and max(abs(residual))<2e-7,**entry)
    # Requested empirical convergence check for the reported OTO curve only.
    oto_convergence=[]
    coarse=Spectrum(Smearing("gaussian",.3),order=20,periods=50)
    fine=Spectrum(Smearing("gaussian",.3),order=40,periods=50)
    oto_clock=matched_clocks(10,.5)["compact"]
    saved_curve=json.loads((OUT/"otoc-curve.json").read_text())
    for t in [.7,1.5,3.]:
        c20=coarse.c4((t,0,t,0),oto_clock)
        c40=fine.c4((t,0,t,0),oto_clock)
        row=min(saved_curve,key=lambda r:abs(r["t"]-t))
        saved=complex(row["clock_real"],row["clock_imag"])
        change=abs(c40-c20)/3
        entry={"t":t,"C4_order20":[c20.real,c20.imag],"C4_order40":[c40.real,c40.imag],
               "normalized_change":change,"saved_curve_difference":abs(saved-c20)/3,
               "scope":"Empirical order convergence at these three times only, not an all-time certificate."}
        oto_convergence.append(entry)
        record("OTO-order20-to40",change<2e-6 and abs(saved-c20)/3<2e-8,**entry)
    finite_mode_check()
    tail_checks=time_tail_checks(cases)
    failures=[r for r in RECORDS if not r["passed"]]
    hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in
            [Path(__file__),ROOT/"code/field_model.py",ROOT/"code/clock_model.py",
             ROOT/"code/run_field_study.py"]}
    output={"status":"passed" if not failures else "failed","checks":len(RECORDS),
            "failed":len(failures),"code_sha256":hashes,"cases_compared":len(comparisons),
            "input_sha256":{p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest()
                            for p in [OUT/"study.json",OUT/"otoc-curve.json"]},
            "runtime_seconds":time.perf_counter()-start,
            "versions":{"python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__},
            "maximum_C2_absolute_difference":max(abs(r["C2_signed_difference"]) for r in comparisons),
            "maximum_C4_absolute_difference":max(abs(r["C4_signed_difference"]) for r in comparisons),
            "comparisons":comparisons,"tail_tasks":tails,"clock_time_tails":tail_checks,
            "CDF_reference_crosschecks":cdf_crosschecks,"OTO_convergence":oto_convergence,
            "failures":failures,"records":RECORDS,"warnings":WARNINGS,
            "coverage":["All 72 cases versus adaptive one-dimensional q convolutions/CDF, not a tensor grid",
                        "Physical spectrum prefactor via Euclidean Wightman function",
                        "Time-smearing normalization/RMS, direct Fourier transform, spectral/time Wightman convolution",
                        "KMS spectral relation and Fourier boundary",
                        "48-level thermal oscillator GNS projected products versus Wick measures",
                        "Representative adaptive CDF tails versus repaired production CDF tails; initial indicator bias saved separately",
                        "Clock timing tails beyond integration endpoint 200"],
            "limitations":["Adaptive quadrature errors are numerical estimates, not interval certificates.",
                           "Results below 1e-12 are subject to absolute, not relative, accuracy.",
                           "Finite spectral truncation is separately bounded; no finite-G or recoil calculation.",
                           "Shared-clock free-field observables only; no claim of interacting de Sitter duality."]}
    OUT.mkdir(exist_ok=True)
    (OUT/"validation.json").write_text(json.dumps(output,indent=2)+"\n",encoding="utf-8")
    (OUT/"tails.json").write_text(json.dumps({"tasks":tails,"source_validation":"validation.json"},indent=2)+"\n",encoding="utf-8")
    print(json.dumps({k:output[k] for k in ["status","checks","failed","cases_compared",
          "maximum_C2_absolute_difference","maximum_C4_absolute_difference","runtime_seconds"]},indent=2),flush=True)
    if failures:
        print(json.dumps([{"name":r["name"],"field_family":r.get("field_family"),
                          "sigma":r.get("field_sigma"),"clock":r.get("clock_family"),
                          "mean":r.get("clock_mean"),"difference":r.get("C4_signed_difference",
                          r.get("C2_signed_difference",r.get("absolute_error")))} for r in failures],indent=2))
    raise SystemExit(1 if failures else 0)

if __name__=="__main__":
    main()

