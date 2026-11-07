"""Fresh all-rank neutral-slope checks; no old driver is executed."""
from pathlib import Path
from fractions import Fraction as F
from math import comb
from datetime import datetime, timezone
import hashlib, json, platform
import mpmath as mp

ROOT = Path(__file__).resolve().parent.parent
DIR = ROOT / 'results/susskind-slater-neutral-slope'
OUT = DIR / 'rank-checks.json'
if OUT.exists():
    raise FileExistsError(OUT)

def choose(n, k):
    return comb(n, k) if 0 <= k <= n else 0

def exact(N, p):
    m = N // 2
    D = choose(m, p) * choose(N-m+p, p)
    omega = [F((choose(N,r)**2-choose(N,r-1)**2)*choose(m-r,m-p)*
               choose(N-p-r,m-p), choose(N,m)*D) for r in range(p+1)]
    mu = [F(choose(m-r,p)*choose(N-m+p-r,p),D) for r in range(p+1)]
    assert sum(omega) == 1 and min(omega) > 0
    avg = sum(w*x for w,x in zip(omega,mu))
    assert sum(w*(x-avg) for w,x in zip(omega,mu)) == 0
    checks = 0
    for d in range(p):
        r = p-d
        a,b = F(r,N-r+1), F(r-1,N-r+2)
        ratio = a*a*(1-b*b)/(1-a*a)*F(m-r+1,d+1)*F(N-p-r+1,N-m-r+1)
        assert omega[r-1]/omega[r] == ratio
        checks += 1
    if p <= N//8:
        for d in range(p+1):
            assert 0 <= p*(mu[p-d]-mu[p]) <= 6*F(p*p,N)*d
    return dict(N=N,p=p,rank_ratio_checks=checks,regression=str(avg),
                mean_deficit=str(sum(w*(p-r) for r,w in enumerate(omega))))

def scaled(ll, p, precision):
    mp.mp.dps = precision
    lam = mp.mpf(ll)
    N = int(p*p/lam)
    assert N%2 == 0 and mp.mpf(p*p)/N == lam
    m=N//2; q=mp.exp(-4*lam); weights=[mp.mpf(1)]
    for d in range(p):
        r=p-d; a=mp.mpf(r)/(N-r+1); b=mp.mpf(r-1)/(N-r+2)
        ratio=a*a*(1-b*b)/(1-a*a)*(mp.mpf(m-r+1)/(d+1))*((mp.mpf(N-p-r+1))/(N-m-r+1))
        assert ratio <= 8*lam/(d+1)
        weights.append(weights[-1]*ratio)
    norm=mp.fsum(weights); weights=[w/norm for w in weights]
    mus=[mp.fprod((1-mp.mpf(p-d)/(m-j))*(1-mp.mpf(p-d)/(N-m+p-j))
                  for j in range(p)) for d in range(p+1)]
    avg=mp.fsum(w*x for w,x in zip(weights,mus))
    alpha=1-mp.mpf(p)/m; eps=1-alpha**2
    pref=m*eps*alpha/(8*lam*q*mp.sqrt(lam))
    fs=[pref*(x-avg) for x in mus]
    gs=[(d-lam)/mp.sqrt(lam) for d in range(p+1)]
    error=mp.fsum(w*(f-g)**2 for w,f,g in zip(weights,fs,gs))
    var=mp.fsum(w*f*f for w,f in zip(weights,fs))
    regression=mp.fsum(w*f for w,f in zip(weights,fs))
    assert abs(regression) < mp.mpf('1e-65')
    # Selected disjoint p-hop has rank p, so its mixed coefficient is exact.
    L=mp.binomial(m,p)**2; D=mp.binomial(m,p)*mp.binomial(m+p,p)
    gammaN=mp.sqrt(L/D)
    mixed=gammaN*fs[0]; target=-mp.sqrt(lam)*mp.exp(-lam)
    return dict(N=N,p=p,lambda_value=ll,precision=precision,
        mean_deficit=mp.nstr(mp.fsum(d*w for d,w in enumerate(weights)),65),
        squared_rank_bridge_error=mp.nstr(error,65), normalized_slope_variance=mp.nstr(var,65),
        regression_residual=mp.nstr(regression,25),
        selected_mixed_coefficient=mp.nstr(mixed,65), mixed_limit=mp.nstr(target,65),
        relative_mixed_error=mp.nstr((mixed-target)/target,40),
        rank_tail_d_ge_12=mp.nstr(mp.fsum(weights[12:]),45),
        all_rank_probabilities=[mp.nstr(w,55) for w in weights],
        all_normalized_slope_coefficients=[mp.nstr(f,55) for f in fs])

finite=[exact(N,p) for N in (8,12,16,24,32) for p in range(1,min(4,N//2)+1)]
rows=[]; changes=[]
for ll in ('0.25','0.5','1.0'):
    for p in (16,32,64,128):
        low=scaled(ll,p,80); high=scaled(ll,p,120)
        change=max(abs(mp.mpf(low[key])-mp.mpf(high[key])) for key in
                   ('squared_rank_bridge_error','normalized_slope_variance','selected_mixed_coefficient'))
        assert change < mp.mpf('1e-60')
        rows.append(high);changes.append(change)
        print(json.dumps({k:high[k] for k in ('p','lambda_value','squared_rank_bridge_error')}),flush=True)
DIR.mkdir(parents=True,exist_ok=True)
inputs=[Path(__file__),ROOT/'results/susskind-slater-covariant/local-evidence.json']
result=dict(status='passed',created_at_utc=datetime.now(timezone.utc).isoformat(),
    finite_exact=finite,rows=rows, max_precision_change=mp.nstr(max(changes),30),
    parameters=dict(lambda_values=['0.25','0.5','1.0'],p_values=[16,32,64,128],
                    precisions=[80,120],normalization='T/(8 lambda_N C_m^(3/2) q_N sqrt(lambda_N))',
                    rank_sum='all p+1 ranks, no cutoff'),
    environment=dict(python=platform.python_version(),mpmath=mp.__version__),
    input_sha256={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
    scope='Exact rank weights/regression/ratio and physical full-rank numerical slope bridge. Analytic dominated-convergence proof is separate. Finite errors retained; no exact-pinching or interval certificate.')
with OUT.open('x',encoding='utf-8') as stream:
    json.dump(result,stream,indent=2);stream.write('\n')
print(json.dumps(dict(status='passed',finite_systems=len(finite),DS_rows=len(rows),sha256=hashlib.sha256(OUT.read_bytes()).hexdigest())))
