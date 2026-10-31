"""Validate modes, independent local integrals, boundary periods and CG matching."""
from pathlib import Path
from datetime import datetime,timezone
from itertools import product
import hashlib,json,platform
import mpmath as mp
import sympy
from ds3_local_vertex import *
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"results/stage4/ds3-local-vertex"

def encode(x):
    if isinstance(x,dict):return {k:encode(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [encode(v) for v in x]
    if isinstance(x,(mp.mpf,mp.mpc)):return mp.nstr(x,65)
    return x

def run():
    checks=[];mode_rows=[];vertex_rows=[];ward_rows=[];direct_rows=[];boundary_rows=[]
    def check(name,observed,expected,tol="1e-45",scale=None):
        scale=max(abs(expected),mp.mpf("1e-70")) if scale is None else mp.mpf(scale)
        error=abs(observed-expected)/scale
        checks.append(dict(name=name,observed=observed,expected=expected,scale=scale,
            scaled_error=error,tolerance=tol,passed=bool(error<mp.mpf(tol))))
    mp.mp.dps=60
    for ns in [".4","1.1","2.3"]:
        nu=mp.mpf(ns)
        for l in [0,1,2,4]:
            for ts in ["-.8","0",".7"]:
                t=mp.mpf(ts);u=mode(nu,l,t)
                du=mp.diff(lambda x:mode(nu,l,x),t)
                ddu=mp.diff(lambda x:mode(nu,l,x),t,2)
                ode=ddu+2*mp.tanh(t)*du+(1+nu**2+l*(l+1)/mp.cosh(t)**2)*u
                kg=1j*mp.cosh(t)**2*(mp.conj(u)*du-mp.conj(du)*u)
                leg=legendre_mode(nu,l,t)
                check(f"KG_{ns}_{l}_{ts}",kg,1)
                check(f"ODE_{ns}_{l}_{ts}",ode,0,scale=1)
                check(f"Legendre_{ns}_{l}_{ts}",u,leg)
                mode_rows.append(dict(nu=nu,l=l,t=t,u=u,KG=kg,ODE=ode,legendre=leg))
    print("modes complete",flush=True)
    samples=[[".7","1.1",".6"],["1.4",".9","1.2"],["3.2","2.7",".5"]]
    patterns=list(product([False,True],repeat=3))
    for strings in samples:
        nus=list(map(mp.mpf,strings))
        for cc in patterns:
            v,cond=vertex(nus,conjugates=cc)
            closed=spherical_vertex(nus,cc);c=reduced_coefficient(nus,cc)
            check(f"000_{strings}_{cc}",v,closed)
            check(f"reduced_norm_{strings}_{cc}",abs(c)**2,reduced_norm(nus,cc))
            vertex_rows.append(dict(nus=nus,conjugates=cc,value=v,closed=closed,
                coefficient=c,norm=abs(c)**2,closed_norm=reduced_norm(nus,cc),cancellation=cond))
        for cc in [(False,False,False),(True,False,False),(True,False,True)]:
            c=reduced_coefficient(nus,cc)
            for l in [1,2,3]:
                v,cond=vertex(nus,(l,l,0),cc);b=boundary_ll0(nus,l,cc)
                check(f"CG_ll0_{strings}_{cc}_{l}",v,c*b,"1e-42")
                boundary_rows.append(dict(nus=nus,conjugates=cc,l=l,bulk=v,
                                         boundary=b,coefficient=c,cancellation=cond))
            for ls in [(1,0,0),(1,1,1),(2,2,1)]:
                terms=[]
                for i,l in enumerate(ls):
                    for shift in [-1,1]:
                        if l+shift<0:continue
                        j=l+1 if shift==1 else l
                        b=mp.sqrt(j*j*(j*j+nus[i]**2)/((2*j-1)*(2*j+1)))
                        other=list(ls);other[i]+=shift
                        val,_=vertex(nus,other,cc)
                        terms.append((-1 if cc[i] else 1)*b*val)
                residual=mp.fsum(terms);scale=mp.fsum(abs(t) for t in terms)
                check(f"boost_Ward_{strings}_{cc}_{ls}",residual,0,"1e-40",scale=scale)
                ward_rows.append(dict(nus=nus,conjugates=cc,ls=ls,terms=terms,
                                      residual=residual,scale=scale))
    print("vertices and covariance complete",flush=True)
    # Independent time quadrature: finite interval + explicit analytic tail.
    for strings,ls,cc in [
            (samples[0],(0,0,0),(True,False,True)),
            (samples[0],(2,2,0),(True,False,False)),
            (samples[1],(2,1,1),(True,False,True)),
            (samples[2],(0,0,0),(True,False,True))]:
        nus=list(map(mp.mpf,strings));T=96
        val,condition=vertex(nus,ls,cc)
        q=direct_vertex(nus,ls,cc,T=T,step=3)
        q2=direct_vertex(nus,ls,cc,T=T,step=2)
        bound=tail_bound(nus,ls,cc,T)
        check(f"quadrature_refinement_{strings}_{ls}",q,q2,"1e-35")
        # Quadrature error is empirical; the tail bound is separately analytic.
        check(f"direct_vs_beta_{strings}_{ls}",q2,val,"1e-34")
        direct_rows.append(dict(nus=nus,ls=ls,conjugates=cc,T=T,step=[3,2],
            beta=val,direct=q,refined=q2,analytic_tail_bound=bound,
            tail_relative=bound/abs(val),relative_difference=abs(q2-val)/abs(val),
            cancellation=condition))
        print("direct integral",strings,ls,flush=True)
    # An independent series at a principal-series point; mpmath hyper uses convergence acceleration.
    for a in [[mp.mpf("0"),mp.mpf(".3"),mp.mpf("-.2")],
              boundary_exponents(list(map(mp.mpf,samples[0])),(True,False,True))]:
        A=sum(a)
        hyper=4**A/mp.fprod(1+x for x in a)*mp.hyper(
            [1,mp.mpf("1.5")]+[-x for x in a],[mp.mpf(".5")]+[2+x for x in a],1)
        gam=sphere_period(a)
        check(f"sphere_Dougall_{len(boundary_rows)}",hyper,gam,"1e-40")
        boundary_rows.append(dict(kind="Dougall_series",a=a,series=hyper,gamma=gam))
    # Independent conventions: J(0)=1, a single edge average, and unitary shadow phases.
    check("sphere_period_zero",sphere_period([0,0,0]),1)
    for a in map(mp.mpf,["-.4",".3","1.2"]):
        q=mp.quad(lambda x:(2-2*x)**a,[-1,0,1])/2
        check("single_edge_"+str(a),sphere_period([a,0,0]),q,"1e-34")
    shadow_rows=[]
    for ns in [".4","1.1"]:
        nu=mp.mpf(ns)
        for l in [0,1,2,4]:
            s=mp.rf(1-1j*nu,l)/mp.rf(1+1j*nu,l)
            ratio=future_phase(-nu,l)/future_phase(nu,l)
            check(f"shadow_phase_{ns}_{l}",s,ratio)
            check(f"shadow_unitarity_{ns}_{l}",abs(s),1)
            shadow_rows.append(dict(nu=nu,l=l,eigenvalue=s,phase_ratio=ratio))
    # Repeat nontrivial coefficient at lower precision without overwriting higher data.
    precision_rows=[]
    for row in boundary_rows:
        if row.get("l")!=3:continue
        with mp.workdps(40):
            v,_=vertex(row["nus"],(3,3,0),row["conjugates"])
        check("precision_l3_"+str(len(precision_rows)),v,row["bulk"],"1e-30")
        precision_rows.append(dict(low_dps=40,high_dps=60,low=v,high=row["bulk"]))
    files=["code/ds3_local_vertex.py","code/validate_ds3_local_vertex.py",
        "sources/lagogiannis-maloney-wang-1106.2846v1.html",
        "sources/marolf-morrison-0810.5163v1.html","sources/dlmf-16.4.html",
        "sources/derkachov-spiridonov-1711.07073v2.html",
        "sources/bros-et-al-0901.4223v1.html","sources/bros-et-al-0901.4223v1.pdf",
        "sources/kolchmeyer-liu-2411.08090v1.html"]
    data=dict(status="passed" if all(r["passed"] for r in checks) else "failed",
        created_at_utc=datetime.now(timezone.utc).isoformat(),python=platform.python_version(),
        mpmath=mp.__version__,sympy=sympy.__version__,working_dps=60,
        checks=checks,modes=mode_rows,vertices=vertex_rows,boundary=boundary_rows,
        ward=ward_rows,direct=direct_rows,shadow=shadow_rows,precision=precision_rows,
        scope="Declared G=0 BD/KG local scalar vertex on finite spherical modes. No physical mass measure, finite-G sewing, or interval quadrature certificate.",
        hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in files})
    OUT.mkdir(parents=True,exist_ok=True);dst=OUT/"diagnostics.json"
    assert not dst.exists(),"Do not overwrite evidence."
    dst.write_text(json.dumps(encode(data),indent=2)+"\n",encoding="utf-8")
    print(json.dumps(dict(status=data["status"],checks=len(checks),
        failed=[r["name"] for r in checks if not r["passed"]])))

if __name__=="__main__":run()
