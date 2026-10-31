"""Independent DS, group-integral, four-body and antipodal-limit diagnostics."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,platform
import mpmath as mp
import numpy as np
from ds3_two_matter import *
from ds3_group_pairing import spherical

ROOT=Path(__file__).resolve().parents[1]

def encode(x):
    if isinstance(x,dict):return {k:encode(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [encode(v) for v in x]
    if isinstance(x,(mp.mpf,mp.mpc)):return mp.nstr(x,60)
    if isinstance(x,complex):return mp.nstr(mp.mpc(x),60)
    return x

def run():
    mp.mp.dps=50
    checks=[];weights=[];products=[];four=[];chains=[];shifted=[];limits=[]
    def check(name,observed,expected,tol="1e-36",scale=None):
        a,b=mp.mpc(observed),mp.mpc(expected)
        scale=max(abs(b),mp.mpf("1e-60")) if scale is None else mp.mpf(scale)
        err=abs(a-b)/scale
        checks.append(dict(name=name,observed=a,expected=b,scale=scale,
                           scaled_error=err,tolerance=tol,passed=bool(err<mp.mpf(tol))))
    for ns in [[".7",".6","1.1"],["1.4","1.2",".9"],["3.2",".5",".4"]]:
        a,q,b=map(mp.mpf,ns)
        group=pair_weight(a,q,b);ds=pair_weight_ds(a,q,b)
        check("DS_pair_"+str(ns),ds,group)
        weights.append(dict(nus=[a,q,b],group=group,ds=ds))
    for ns in [[".7",".6"],["1.4","1.2"]]:
        a,q=map(mp.mpf,ns);L=32
        tail=2*mp.sinh(mp.pi*a)/a*mp.sinh(mp.pi*q)/q*mp.exp(-mp.pi*L)*(L/mp.pi+1/mp.pi**2)
        for r in [mp.mpf(0),mp.mpf(".7"),mp.mpf("2")]:
            value=mp.quad(lambda b:pair_weight(a,q,b)*spherical(b,r),list(range(0,L+1,2)))
            target=spherical(a,r)*spherical(q,r)
            check("spherical_product_"+str(ns)+"_"+str(r),value,target)
            check("spherical_product_tail_"+str(ns)+"_"+str(r),value-target,0,"1",tail+mp.mpf("1e-46"))
            products.append(dict(nus=[a,q],r=r,L=L,direct=value,target=target,tail=tail))
    print("DS normalization and spherical completeness complete",flush=True)
    for ns in [[".7","1.1",".6","1.2"],["1","1","1","1"],["1.4",".9","1.2",".5"]]:
        a,q,c,p=map(mp.mpf,ns);nus=[a,q,c,p];T=40;L=18
        closed=four_spherical_closed(nus)
        radial=mp.quad(lambda r:mp.fprod(mp.sin(n*r) for n in nus)/mp.sinh(r)**2 if r else 0,
                       list(range(0,T+1,2)))/mp.fprod(nus)
        spectral=mp.quad(lambda b:four_weight(a,b,c,p,q),list(range(0,L+1,2)))
        rtail=2/(mp.fprod(nus)*mp.expm1(2*T))
        stail=mp.fprod(mp.sinh(mp.pi*n)/n for n in nus)*mp.exp(-2*mp.pi*L)
        check("four_radial_"+str(ns),radial,closed,"1e-30")
        check("four_spectral_"+str(ns),spectral,closed,"1e-35")
        check("four_radial_tail_"+str(ns),radial-closed,0,"1",rtail+mp.mpf("1e-45"))
        check("four_spectral_tail_"+str(ns),spectral-closed,0,"1",stail+mp.mpf("1e-45"))
        four.append(dict(nus=nus,T=T,L=L,closed=closed,radial=radial,spectral=spectral,
                         radial_tail=rtail,spectral_tail=stail))
    for ns in [[".7","1.1",".9",".6","1.2"],["2.3",".4","1.8",".5",".7"]]:
        a,b,c,p,q=map(mp.mpf,ns)
        amp=mp.sech(mp.pi*c)*emission000(b,c,p)*emission000(a,b,q)
        norm=four_weight(a,b,c,p,q);raw=abs(amp)**2/norm;expected=chain(a,b,c,p,q)
        check("chain_quotient_"+str(ns),raw,expected)
        chains.append(dict(nus=[a,b,c,p,q],amplitude=amp,rigging_weight=norm,raw=raw,expected=expected))
    for ns in [[".4","1.2",".6"],["2","3",".6"],["5","4","1.2"]]:
        a,b,p=map(mp.mpf,ns)
        delta=mp.sqrt(1+a*a)-mp.sqrt(1+b*b)
        raw=joint(a,b,p)/state_density(b)*mp.exp(mp.pi*delta)
        value=tilted(a,b,p)
        check("tilt_identity_"+str(ns),raw,value)
        shifted.append(dict(nus=[a,b,p],energy=delta,raw=raw,tilted=value,g_ratio=clock_g(a)/clock_g(b),
                            balanced=balanced(a,b,p),limiting_spectrum=k(a-b,p)))
    print("four-body norm and exact half-period identities complete",flush=True)
    s=.2;p=.6;q=1.2;D=14
    for M in [4,8,16,32]:
        for tp,tq in [(0,0),(.5,1),(2,.7)]:
            coarse=quadrature(M,s,p,q,tp,tq,64,D)
            refined=quadrature(M,s,p,q,tp,tq,96,D)
            pm,qm=map(mp.mpf,map(str,[p,q]));tm,tn=map(mp.mpf,map(str,[tp,tq]))
            leading=lr(tm,pm)*lr(tn,qm)
            derivative=mp.diff(lambda t:lr(t,pm),tm)*lr(tn,qm)+lr(tm,pm)*mp.diff(lambda t:lr(t,qm),tn)
            first=1j*derivative/M
            second=(-1j*mp.pi*derivative/2-
                    (tm*mp.diff(lambda t:lr(t,pm),tm)*lr(tn,qm)+
                     tn*lr(tm,pm)*mp.diff(lambda t:lr(t,qm),tn))/2)/M**2
            bound=limit_bound(mp.mpf(M),mp.mpf(str(s)),pm,qm,max(abs(tm),abs(tn)))
            tail=tail_bound(mp.mpf(D),mp.mpf(str(s)))
            check(f"four_time_refine_{M}_{tp}_{tq}",coarse,refined,"2e-10")
            check(f"four_limit_bound_{M}_{tp}_{tq}",mp.mpc(refined)-leading,0,"1",bound)
            limits.append(dict(M=M,cutoff=s,p=p,q=q,tp=tp,tq=tq,D=D,orders=[64,96],
                coarse=coarse,refined=refined,leading=leading,first=first,second_formal=second,
                error=abs(mp.mpc(refined)-leading),relative_error=abs(mp.mpc(refined)/leading-1),
                corrected_error=abs(mp.mpc(refined)-leading-first),uniform_absolute_bound=bound,outer_tail=tail))
        print("proper-time two-matter quadratures M =",M,flush=True)
    paths=["code/ds3_two_matter.py","code/validate_ds3_two_matter.py","code/ds3_group_pairing.py",
        "code/ds3_local_vertex.py","code/ds3_mass_spectrum.py",
        "sources/derkachov-spiridonov-1711.07073v2.html",
        "sources/kolchmeyer-liu-2411.08090v1.html","sources/marolf-morrison-0810.5163v1.html",
        "sources/marolf-yu-2610.24073v1.html"]
    data=dict(status="passed" if all(c["passed"] for c in checks) else "failed",
        created_at_utc=datetime.now(timezone.utc).isoformat(),python=platform.python_version(),
        mpmath=mp.__version__,numpy=np.__version__,dps=50,checks=checks,
        weights=weights,products=products,four_spherical=four,chains=chains,shifted=shifted,limits=limits,
        scope="Two distinct species, ordered scalar observer channel; exact physical norm and Euclidean-shifted uncrossed benchmark. Not crossed recoupling, same-species Fock algebra or finite gravity. No interval certificate; second_formal is diagnostic only.",
        hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths})
    out=ROOT/"results/stage4/ds3-two-matter";out.mkdir(parents=True,exist_ok=True)
    dest=out/"diagnostics.json";assert not dest.exists()
    dest.write_text(json.dumps(encode(data),indent=2)+"\n",encoding="utf-8")
    print(json.dumps(dict(status=data["status"],checks=len(checks),failed=[c["name"] for c in checks if not c["passed"]])))
if __name__=="__main__":run()
