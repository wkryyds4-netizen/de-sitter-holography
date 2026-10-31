"""Finite preparation four-point study, with independent mass-variable reference."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,platform
import numpy as np
import scipy
from scipy.integrate import quad
from scipy.special import loggamma
from quantized_observer import DoublePacket, mass, s_of_p, rho_pl, hh_kernel, right_kernel
from observer_fourpoint import FourPointQuadrature

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"results/stage3/fourpoint"
checks=[]
def serial(x):
    if isinstance(x,np.ndarray): return serial(x.tolist())
    if isinstance(x,np.generic): return serial(x.item())
    if isinstance(x,complex): return {"real":x.real,"imag":x.imag}
    if isinstance(x,dict): return {str(k):serial(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [serial(v) for v in x]
    return x
def check(name,passed,**details):
    checks.append(dict(name=name,passed=bool(passed),**serial(details)))
def close(name,x,y,atol=2e-10,rtol=2e-8):
    e=float(np.max(np.abs(np.asarray(x)-np.asarray(y))))
    scale=float(np.max(np.abs(np.asarray(y))))
    check(name,e<=atol+rtol*scale,error=e,tolerance=atol+rtol*scale)
    return e
def raw_hh(s):
    return float(rho_pl(s))/np.cosh(np.pi*s)**2

def independent_mass_h(m_initial,nu,sigma,m_lower,u=0.):
    """Integrate physical final mass, not spectral label s; a tail to infinity.
    Computes Gamma8 afresh from all eight factors, not right_kernel().
    """
    b=np.sqrt(m_initial*m_initial-.25)
    def integrand(ma):
        a=np.sqrt(max(0.,ma*ma-.25))
        if a==0:return 0j
        lg=sum(loggamma(.25+.5j*(sn*nu+sa*a+sb*b))
               for sn in (-1,1) for sa in (-1,1) for sb in (-1,1))
        logcosh_b=np.logaddexp(np.pi*b,-np.pi*b)-np.log(2)
        log_weight=(np.log(a*np.tanh(np.pi*a))+2*logcosh_b
                    -np.log(16)-5*np.log(np.pi)+lg.real
                    -sigma*sigma*(ma-m_initial)**2+np.log(ma/a))
        return np.exp(log_weight+1j*u*(ma-m_initial))
    breaks=sorted(set([m_lower]+[x for x in (m_initial-nu,m_initial,m_initial+nu,m_initial+8) if x>m_lower]))
    val=0j;err=0.
    for left,right in zip(breaks,breaks[1:]+[np.inf]):
        re,er=quad(lambda z:integrand(z).real,left,right,epsabs=2e-11,epsrel=2e-10)
        im,ei=quad(lambda z:integrand(z).imag,left,right,epsabs=2e-11,epsrel=2e-10)
        val+=re+1j*im;err+=er+ei
    return val,err

def main():
    packet=DoublePacket()
    tgrid=[-2.,-1.,0.,1.,2.]
    tuples=[[t,0.,t,0.] for t in tgrid]
    rows=[]
    max_refine=max_extend=0.
    for m0 in (.5,2.,5.):
        for nu in (.5,1.):
            for sigma in (.25,.5):
                for projected in (True,False):
                    if m0==.5 and not projected:continue # identical lower limit, checked below
                    args=(packet,m0,nu,sigma,projected)
                    coarse=FourPointQuadrature(*args,p_order=24,mass_order=16,tail_extra=10)
                    fine=FourPointQuadrature(*args,p_order=48,mass_order=28,tail_extra=10)
                    extended=FourPointQuadrature(*args,p_order=48,mass_order=28,tail_extra=14)
                    vals=[fine.values(t) for t in tuples]
                    z=fine.zero_diagnostics()
                    tag=f"{m0=},{nu=},{sigma=},{projected=}"
                    for channel in ("AABB","ABBA"):
                        mainvals=[v[channel] for v in vals]
                        err=close("Gauss refinement "+channel+" "+tag,
                                  [coarse.values(t)[channel] for t in tuples],mainvals)
                        ext=close("cutoff extension "+channel+" "+tag,
                                  [extended.values(t)[channel] for t in tuples],mainvals)
                        max_refine=max(max_refine,err);max_extend=max(max_extend,ext)
                        close("real-time conjugation "+channel+" "+tag,mainvals,np.conj(mainvals[::-1]))
                        check("zero-time positive "+channel+" "+tag,z[channel]>0,value=z[channel])
                        check("all sampled time magnitudes bounded by zero "+channel+" "+tag,
                              max(abs(x) for x in mainvals)<=z[channel]+2e-10)
                        close("common time translation "+channel+" "+tag,
                              fine.values([x+.37 for x in tuples[-1]])[channel],mainvals[-1])
                    close("adjacent variance identity "+tag,z["adjacent_excess"],z["variance_of_conditional_C2"])
                    check("adjacent Jensen bound "+tag,z["adjacent_excess"]>=-2e-11,
                          value=z["adjacent_excess"])
                    rows.append(dict(parameters=fine.parameters,times_left_to_right=tuples,values=vals,
                                     zero=z,model="two independent species of equal mass",
                                     scope="each AABB/ABBA is a complete mixed-species moment; ABAB not computed here"))
                    print(tag,"adjacent/nested normalized",z["normalized_adjacent"],z["normalized_nested"],flush=True)
    for m0 in (2.,5.):
        for nu in (.5,1.):
            for sigma in (.25,.5):
                select=[r for r in rows if (r["parameters"]["m0"],r["parameters"]["nu"],r["parameters"]["sigma_proper"])==(m0,nu,sigma)]
                projected=next(r for r in select if r["parameters"]["projected"])
                unrestricted=next(r for r in select if not r["parameters"]["projected"])
                for k in ("AABB","ABBA"):
                    check(f"zero-time cutoff monotonicity {m0,nu,sigma,k}",
                          projected["zero"][k]<=unrestricted["zero"][k]+2e-10)
    # Independently verify HH diagram weights to general packet conversion.
    for d,a,b,c in [(1.1,.8,1.9,2.1),(3.2,2.8,4.1,3.5)]:
        for nu in (.5,1.):
            raw_nested=hh_kernel(nu,a,b)*hh_kernel(nu,a,d)/raw_hh(a)/raw_hh(d)
            factored_nested=right_kernel(nu,a,d)*right_kernel(nu,b,a)
            close(f"6.23 HH conversion {nu,d,a,b}",raw_nested,factored_nested,atol=1e-15,rtol=1e-11)
            raw_adjacent=hh_kernel(nu,d,c)*hh_kernel(nu,d,a)/raw_hh(d)**2
            close(f"6.24 HH conversion {nu,d,a,c}",raw_adjacent,
                  right_kernel(nu,a,d)*right_kernel(nu,c,d),atol=1e-15,rtol=1e-11)
    # Actual independent reference in mass variables, at two fixed exterior masses.
    reference=[]
    refq=FourPointQuadrature(packet,2.,1.,.5,True,p_order=48,mass_order=28,tail_extra=14)
    for idx in (12,70):
        b=float(refq.b[idx]);mb=float(mass(b))
        h,err=independent_mass_h(mb,1.,.5,2.,.7)
        close(f"independent eight-Gamma mass-integral h at b={b}",refq.h_external(.7)[idx],h,
              atol=2e-9,rtol=2e-8)
        reference.append(dict(b=b,m=mb,u=.7,value=h,quad_error_estimate=err,
                              production=refq.h_external(.7)[idx]))
    # Full adjacent packet integral, independently adaptive in p and final physical mass.
    print("Independent adjacent four-point mass-variable integration",flush=True)
    max_inner=[0.]
    def adj_integrand(p):
        mb=2.+p/(2*np.pi)
        val,err=independent_mass_h(mb,1.,.5,2.)
        max_inner[0]=max(max_inner[0],err)
        return float(packet.probability(p))*val.real**2
    ref_adj=0.;err_adj=0.
    for lo,hi in packet.intervals:
        value,error=quad(adj_integrand,lo,hi,epsabs=2e-10,epsrel=2e-8)
        ref_adj+=value;err_adj+=error
    close("independent full adjacent four-point",refq.aabb([0,0,0,0]),ref_adj,atol=2e-9,rtol=2e-8)
    # Nested reference: independent physical mass integrals for fixed exterior b.
    # The outer packet quadrature is separately checked above; do not call this a full adaptive packet reference.
    print("Independent nested conditional mass-variable integration",flush=True)
    fixed_mb=2.+6/(2*np.pi)
    fixed_b=float(np.sqrt(fixed_mb**2-.25))
    def nested_outer(ma):
        a=np.sqrt(max(0.,ma*ma-.25))
        if a==0:return 0.
        h,_=independent_mass_h(ma,1.,.5,2.)
        return float(right_kernel(1.,a,fixed_b))*np.exp(-.25*(ma-fixed_mb)**2)*(ma/a)*h.real
    ref_n,err_n=quad(nested_outer,2.,fixed_mb+14,points=[fixed_mb],epsabs=2e-9,epsrel=2e-7,limit=100)
    delta=mass(refq.x)-fixed_mb
    # Get quadrature weights without accessing internal grid transition.
    from quantized_observer import composite_gauss
    xx,ww=composite_gauss(refq.lower,refq.upper,28)
    cond=float(np.sum(ww*right_kernel(1.,xx,fixed_b)*np.exp(-.25*delta**2)*refq.h_grid(0).real))
    close("independent nested conditional integral",cond,ref_n,atol=2e-8,rtol=2e-7)
    scripts=["observer_fourpoint.py","run_observer_fourpoint.py","quantized_observer.py"]
    hashes={name:hashlib.sha256((ROOT/"code"/name).read_bytes()).hexdigest() for name in scripts}
    meta=dict(generated_utc=datetime.now(timezone.utc).isoformat(),code_sha256=hashes,
              versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
              source="https://arxiv.org/html/2411.08090v1",source_equations=["5.13","6.10","6.23","6.24","6.25"],
              units="L=1; sigma and four times proper time; p=2*pi*(m-m0)",
              limitations=["G->0, dS2, two independent same-mass real principal-series species",
                           "Not a same-species full four-point or the crossed ABAB OTOC",
                           "Gauss refinement and cutoff differences are empirical; no interval quadrature certification"])
    study=dict(**meta,packet=dict(c1=3,c2=9,width=1,mean_p=6),rows=rows,
               references=dict(conditional_h=reference,adjacent_packet_value=ref_adj,
                               adjacent_outer_error_estimate=err_adj,maximum_inner_error_estimate=max_inner[0],
                               nested_fixed_exterior_mass=fixed_mb,nested_value=ref_n,
                               nested_outer_error_estimate=err_n,nested_production=cond,
                               nested_reference_cutoff=fixed_mb+14,
                               nested_scope="conditional kernel; finite outer cutoff and infinite inner mass integral"))
    failed=sum(not c["passed"] for c in checks)
    validation=dict(**meta,status="passed" if failed==0 else "failed",checks=len(checks),failed=failed,
                    check_details=checks,maxima=dict(refinement=max_refine,cutoff_extension=max_extend))
    OUT.mkdir(parents=True,exist_ok=True)
    for name,data in [("uncrossed-study.json",study),("uncrossed-validation.json",validation)]:
        (OUT/name).write_text(json.dumps(serial(data),ensure_ascii=False,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps({"status":validation["status"],"checks":len(checks),"failed":failed,"maxima":validation["maxima"]}),flush=True)
    if failed:raise SystemExit(1)

if __name__=="__main__":main()
