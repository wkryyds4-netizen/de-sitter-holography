"""Independent finite five-word checks of the central-window filter drift."""
from pathlib import Path
from math import comb
from fractions import Fraction as F
import hashlib,json,platform
import sympy as s
import mpmath as mp
from memory_two_charge_filter import filters

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-central-window'


def main():
    target=OUT/'central-filter-checks.json'
    if target.exists():raise FileExistsError(target)
    old=json.loads((ROOT/'results/susskind-memory-two-charge/symbolic-filter.json').read_text())
    apre=json.loads((ROOT/'results/susskind-memory-filter-limit/filter-certificate.json').read_text())
    hcoeff=list(map(s.sympify,old['H_squared_coefficients']))
    acoeff=list(map(s.sympify,apre['exact_prefactors']))
    m,n,p,a,b,r=s.symbols('m n p a b r')
    records=[];checks=0
    for pp in [16,32,64,128,256]:
        nn=4*pp*pp
        for xx in [-1,0,1]:
            mm=nn//2+xx*pp
            for u,v,S in [(0,0,0),(2,1,3)]:
                aa,bb=pp-u,pp-v;rr=aa+bb-S
                C=comb(mm,pp)*comb(nn-mm+pp,pp)
                vals={m:mm,n:nn-mm+1,p:pp,a:aa,b:bb,r:rr}
                kap=lambda j:comb(mm-j,pp)*comb(nn-mm+pp-j,pp)
                g=filters(nn,pp,mm,aa,bb,rr)['g']
                gh=g[0]+(mm-2*pp)*g[1]
                for coeff,want in [(hcoeff,gh),(acoeff,g[1])]:
                    got=sum(c.subs(vals)*k for c,k in zip(coeff,[C,kap(aa),kap(bb),kap(rr)]))
                    assert got==s.Rational(want.numerator,want.denominator)
                    checks+=1
                mh=nn//2;Ch=comb(mh,pp)*comb(nn-mh+pp,pp)
                ghalf=filters(nn,pp,mh,aa,bb,rr)['g']
                hhalf=ghalf[0]+(mh-2*pp)*ghalf[1]
                shift=gh/C-hhalf/Ch
                runs=[]
                for digits in [80,120]:
                    mp.mp.dps=digits
                    cv=lambda z:mp.mpf(z.numerator)/z.denominator
                    q=mp.exp(-1);pred=xx*(1-q)**2
                    runs.append([cv(shift),pred,cv(shift)-pred])
                assert all(abs(lo-hi)<mp.mpf('1e-70') for lo,hi in zip(*runs))
                got,pred,dev=runs[-1]
                records.append(dict(N=nn,p=pp,x=xx,deficits=[u,v,S],
                    exact_central_minus_half=str(shift),central_minus_half=mp.nstr(got,65),
                    predicted_limit=mp.nstr(pred,65),absolute_limit_deviation=mp.nstr(dev,45),
                    relative_limit_deviation=mp.nstr(dev/pred,45) if xx else None))
    names=['code/check_memory_central_filter.py','code/memory_two_charge_filter.py',
           'results/susskind-memory-two-charge/symbolic-filter.json',
           'results/susskind-memory-filter-limit/filter-certificate.json']
    data=dict(status='passed',exact_filter_identities=checks,records=records,precisions=[80,120],
        scope='Exact finite rational agreement with independent five-word coefficients. Convergence samples do not replace the analytic asymptotic expansion; zero-limit cases use absolute error.',
        environment=dict(python=platform.python_version(),sympy=s.__version__,mpmath=mp.__version__),
        input_sha256={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in names})
    target.write_text(json.dumps(data,indent=2)+'\n')
    print('passed',checks,'exact identities; last finite shift',records[-1]['central_minus_half'],flush=True)


if __name__=='__main__':main()
