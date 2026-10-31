"""Non-spherical four-body norm: pair l=1 singlet against four spherical legs."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,platform
import mpmath as mp
from ds3_local_vertex import vertex
from ds3_group_pairing import spherical
from ds3_two_matter import four_weight
from validate_ds3_group_offdiagonal import boost_l0
ROOT=Path(__file__).resolve().parents[1]

def run():
    mp.mp.dps=50;checks=[];rows=[];vertices=[]
    def check(name,a,b,tol="1e-35",scale=None):
        scale=max(abs(b),mp.mpf("1e-60")) if scale is None else scale
        error=abs(a-b)/scale
        checks.append(dict(name=name,observed=mp.nstr(a,60),expected=mp.nstr(b,60),
            scale=mp.nstr(scale,60),scaled_error=mp.nstr(error,60),tolerance=tol,passed=bool(error<mp.mpf(tol))))
    for ss in [[".7","1.1",".6","1.2"],["1.4",".9","1.2",".5"]]:
        a,q,c,p=map(mp.mpf,ss);T=40;L=18
        A=1+a*a+q*q;den=2*mp.sqrt((1+a*a)*(1+q*q))
        def inv(b):return mp.sqrt(3)*(A-b*b)/den
        for b in map(mp.mpf,[".4","1","2"]):
            v0=vertex([a,q,b],(0,0,0),(True,False,True))[0]
            v1=vertex([a,q,b],(1,1,0),(True,False,True))[0]
            value=mp.sqrt(3)*v1/v0
            check("l1_invariant_"+str(ss)+"_"+str(b),value,inv(b))
            vertices.append(dict(nus=ss,b=str(b),vertex=mp.nstr(value,60),polynomial=mp.nstr(inv(b),60)))
        direct=mp.quad(lambda r:-mp.sinh(r)**2/mp.sqrt(3)*boost_l0(a,1,r)*boost_l0(q,1,r)*
            spherical(c,r)*spherical(p,r),list(range(0,T+1,2)))
        spectral=mp.quad(lambda b:four_weight(a,b,c,p,q)*inv(b),list(range(0,L+1,2)))
        rtail=2*mp.sqrt(3)/(c*p*(1-mp.exp(-2*T))**2)*mp.exp(-2*T)*(T*T+T+mp.mpf(".5"))
        beta=2*mp.pi
        stail=2*mp.pi*mp.sqrt(3)/den*mp.fprod(mp.sinh(mp.pi*n)/n for n in [a,q,c,p])*mp.exp(-beta*L)*(
            (L*L+A)/beta+2*L/beta**2+2/beta**3)
        check("four_l1_"+str(ss),direct,spectral,"1e-29")
        check("four_l1_tail_"+str(ss),direct-spectral,0,"1",rtail+stail+mp.mpf("1e-43"))
        rows.append(dict(nus=ss,T=T,L=L,direct=mp.nstr(direct,60),spectral=mp.nstr(spectral,60),
                         radial_tail=mp.nstr(rtail,60),spectral_tail=mp.nstr(stail,60)))
    paths=["code/validate_ds3_two_matter_angular.py","code/ds3_two_matter.py",
           "code/validate_ds3_group_offdiagonal.py","code/ds3_local_vertex.py","code/ds3_group_pairing.py"]
    data=dict(status="passed" if all(c["passed"] for c in checks) else "failed",
        created_at_utc=datetime.now(timezone.utc).isoformat(),python=platform.python_version(),
        mpmath=mp.__version__,dps=50,checks=checks,rows=rows,vertices=vertices,
        scope="Two independent non-spherical four-body group integrals and six l=1 invariant checks; analytic outer tails, no interval certification; not a spin-summed crossing test.",
        hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths})
    dest=ROOT/"results/stage4/ds3-two-matter/angular.json";assert not dest.exists()
    dest.write_text(json.dumps(data,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(dict(status=data["status"],checks=len(checks),failed=[c["name"] for c in checks if not c["passed"]])))
if __name__=="__main__":run()
