"""Full-rank neutral moment and candidate coherent cross; no parity rule."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,platform,time
import mpmath as mp
from cubic_linear_limit import distribution
from cubic_linear_dynamics import weights,crossing
from memory_central_norm import exchanged
from memory_central_cross import pair_cross
from memory_coherent_cross import charged

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-central-cross'

def evaluate(N,p,m,dps):
    mp.mp.dps=dps
    w=list(reversed(distribution(N,p,m)))
    mu=mp.fsum(w[a]*(p-a) for a in range(p+1))
    def avg(z,centered=False):
        return mp.fsum(w[a]*z[a]*((p-a-mu) if centered else 1) for a in range(p+1))
    aa=[mp.mpf(a) for a in range(p+1)];xx=[mp.mpf(a*(N+1-a)) for a in range(p+1)]
    ma=avg(aa);mx=avg(xx);mac=avg(aa,True);mxc=avg(xx,True)
    sums=[[] for _ in range(6)];zz=[mp.mpf(1)]*(p+1);hh=mp.mpf(1)
    for j in range(p+1):
        z0=avg(zz);za=avg([a*z for a,z in zip(aa,zz)])
        zac=avg([a*z for a,z in zip(aa,zz)],True)
        hd=hh*mp.mpf(j)/(m-j)
        for row,val in zip(sums,[hh*z0*z0,hh*za*za,hd*z0*z0,hd*za*za,2*hh*za*zac,2*hd*za*zac]):row.append(val)
        if j<p:
            zz=[z*(-a+j)*(a-N-1+j) for a,z in enumerate(zz)]
            hh/=(-m+j)*(m-N+j)*(-N-1+j)*(j+1)
    chi,abchi,D,abD,abchiL,abDL=map(mp.fsum,sums)
    cf=mp.mpf(m*(N-m)*(N+1))/N;g=mp.mpf((m-p)**2)/m;eta=mp.mpf(p*p)/m
    neutral=-2*(ma*ma+abchi)+mx*mx/cf
    XB=neutral/N
    XL=(-4*ma*mac-2*abchiL+2*mx*mxc/cf)/N
    parts=[eta*neutral/N,g*mx*mx/(N*m*(m-1)*(N+1)),g*(m*D-2*abD)/N]
    XO=mp.fsum(parts);lam=mp.mpf(p*p)/N;q=mp.exp(-4*lam)
    targets=[2*lam*(1-q),mp.mpf(0),4*lam*lam*(1+q)+2*lam*(1-q)]
    vals=[XB,XL,XO]
    return dict(N=N,p=p,m=m,x=str(F(2*m-N,2*p)),precision=dps,
        values=[mp.nstr(v,65) for v in vals],targets=[mp.nstr(v,65) for v in targets],
        signed_limit_deviations=[mp.nstr(v-t,40) for v,t in zip(vals,targets)],
        candidate_cross_cancellation_factor=mp.nstr(mp.fsum(map(abs,parts))/abs(XO),30))

def main():
    target=OUT/'cross-limit.json'
    if target.exists():raise FileExistsError(target)
    start=time.time();identities=0;finite=[]
    for N in range(6,14):
        for m in range(2,N-1):
            for a in range(min(3,m-1,N-m)+1):
                for b in range(min(3,m-1,N-m)+1):
                    D=exchanged(N,m,a,b)-crossing(N,m,a,b)
                    val=F(b*(N+1-b),N+1-2*b)*(charged(N,m,a,b-1)-charged(N,m,a,b)) if b else F(0)
                    assert val==-m*D;identities+=1
    for N,p,m in [(5,1,2),(5,2,3),(6,1,4),(6,2,4)]:
        w=weights(N,m,p);mu=sum(w[a]*(p-a) for a in range(p+1))
        exact=[F(0)]*3
        for a in range(p+1):
            for b in range(p+1):
                xb,xo=pair_cross(N,p,m,a,b);ww=w[a]*w[b]
                exact[0]+=ww*xb;exact[1]+=ww*(2*p-a-b-2*mu)*xb;exact[2]+=ww*xo
        got=evaluate(N,p,m,120)
        for val,ex in zip(got['values'],exact):assert abs(mp.mpf(val)-mp.mpf(ex.numerator)/ex.denominator)<mp.mpf('1e-60')
        finite.append(dict(N=N,p=p,m=m,exact_values=list(map(str,exact))))
    rows=[];changes=[]
    for lam in [F(1,10),F(1,4),F(1,2)]:
        for p in [16,32,64,128]:
            N=int(p*p/lam)
            for x in [-2,0,2]:
                m=N//2+x*p;lo=evaluate(N,p,m,80);hi=evaluate(N,p,m,120)
                change=max(abs(mp.mpf(a)-mp.mpf(b)) for a,b in zip(lo['values'],hi['values']))
                assert change<mp.mpf('1e-55');changes.append(change);rows.append(hi)
            print('central cross full rank',str(lam),p,'passed',flush=True)
    names=['code/check_memory_central_cross_limit.py','code/memory_central_cross.py','code/memory_central_norm.py','code/memory_coherent_cross.py','code/cubic_linear_limit.py','code/cubic_linear_dynamics.py']
    data=dict(status='passed',labels=['X/N--B','X/N--L_input','candidate X/N--O'],contiguous_identities=identities,
        finite_checks=finite,records=rows,max_precision_change=mp.nstr(max(changes),30),
        scope='Full p+1 ranks and terminating sums, without numerical tail truncation. Exact contiguous identities and finite sum matches. The third column remains conditional on the general-charge exchanged Casimir recoupling lemma; numerical checks do not prove that lemma. No interval rounding certification.',
        environment=dict(python=platform.python_version(),mpmath=mp.__version__),elapsed_seconds=time.time()-start,
        input_sha256={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in names})
    target.write_text(json.dumps(data,indent=2)+'\n')

if __name__=='__main__':main()
