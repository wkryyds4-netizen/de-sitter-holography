"""Full-rank growing-window checks and binomial outer-tail suppression."""
from pathlib import Path
import hashlib,json,platform,math
import mpmath as mp
from check_memory_central_norm_limit import evaluate
from check_memory_central_raw import calculate

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-charge-sum'

def tail(p,dps):
    with mp.workdps(dps):
        N=4*p*p;A=mp.mpf(80);R=mp.sqrt(A*mp.log(p));k=int(mp.floor(N/2+p*R))+1
        if k>N:return dict(p=p,N=N,R=str(R),first_upper_charge=k,tail='0',bound='0',N8_tail='0')
        term=mp.exp(mp.loggamma(N+1)-mp.loggamma(k+1)-mp.loggamma(N-k+1)-N*mp.log(2))
        terms=[term]
        for j in range(k,N):
            term*=mp.mpf(N-j)/(j+1);terms.append(term)
        prob=2*mp.fsum(terms);bound=2*mp.exp(-2*(mp.mpf(k)-N/2)**2/N)
        assert prob<=bound and N**8*bound<=2*N**8*p**(-40)
        return dict(p=p,N=N,R=mp.nstr(R,65),first_upper_charge=k,tail=mp.nstr(prob,65),
                    bound=mp.nstr(bound,65),N8_tail=mp.nstr(N**8*prob,65))

def main():
    target=OUT/'growing-window-checks.json'
    if target.exists():raise FileExistsError(target)
    rows=[];changes=[]
    for lam in ['0.25','0.5']:
        for p in [32,64,128]:
            N=int(p*p/mp.mpf(lam));j=math.ceil(math.sqrt(4*math.log(p)))
            for x in [-j,j]:
                lo=evaluate(N,p,N//2+x*p,80);hi=evaluate(N,p,N//2+x*p,120)
                fl=calculate(lam,p,x,80);fh=calculate(lam,p,x,120)
                with mp.workdps(130):
                    errs=[abs(mp.mpf(lo[k])-mp.mpf(hi[k])) for k in ['weighted_A_norm','A_relative_limit_error']]
                    errs+=[abs(mp.mpf(fl[k])-mp.mpf(fh[k])) for k in ['variance_over_C3','normalized_F_proxy_L2_error_squared']]
                    change=max(errs);assert change<mp.mpf('1e-60');changes.append(change)
                rows.append(dict(norm=hi,linear_filter=fh,window='x=+/-ceil(sqrt(4 log p))'))
            print('growing full-rank',lam,p,'passed',flush=True)
    tails=[]
    for p in [16,32,64,128]:
        lo=tail(p,80);hi=tail(p,120)
        with mp.workdps(130):assert abs(mp.mpf(lo['tail'])/mp.mpf(hi['tail'])-1)<mp.mpf('1e-60')
        tails.append(hi)
    names=['code/memory_charge_growing_window.py','code/check_memory_central_norm_limit.py',
           'code/check_memory_central_raw.py','code/cubic_linear_limit.py']
    target.write_text(json.dumps(dict(status='passed',records=rows,binomial_tail=tails,
        precisions=[80,120],saved_digits=65,max_absolute_saved_precision_change=mp.nstr(max(changes),25),
        scope='No rank cutoff, moderate-window norm diagnostics and positive full binomial tail sums. Numerical precision checks are not interval certificates or the proof of uniformity.',
        environment=dict(python=platform.python_version(),mpmath=mp.__version__),
        input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names}),indent=2)+'\n')
    print('saved 12 full-rank groups and four outer tails')

if __name__=='__main__':main()
