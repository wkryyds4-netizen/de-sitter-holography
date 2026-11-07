"""Independent Fourier quadrature of exact finite-support charge coefficients."""
from pathlib import Path
from datetime import datetime,timezone
from math import comb,prod
import hashlib,json,platform,time
import mpmath as mp
from central_ordinary_sources import central_singleton_coefficient

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-central-ordinary/fourier-checks.json'

def calculate(p,internals,precision):
    N=4*p*p;k=max(1,sum(internals));L=N-2*k*p;r=len(internals);T=16
    exact=central_singleton_coefficient(N,p,k,internals)
    with mp.workdps(precision):
        den=2*mp.pi*p*mp.mpf(comb(L,L//2))/mp.power(2,L)
        def integrand(t):
            return mp.cos(t/(2*p))**L*mp.re((p*mp.expm1(1j*t/p))**r)*prod(internals)
        value=2*mp.quad(integrand,[0,1,2,4,8,T])/den
        absolute=2*mp.quad(lambda t:abs(integrand(t)),[0,1,2,4,8,T])/den
        a=mp.mpf(L)/(8*p*p)
        tail=prod(internals)*a**(-mp.mpf(r+1)/2)*mp.gammainc(mp.mpf(r+1)/2,a*T*T,mp.inf)/den
        target=mp.mpf(exact.numerator)/exact.denominator
        return {'value':mp.nstr(value,precision),'analytic_cutoff_tail':mp.nstr(tail,precision),
                'absolute_error':mp.nstr(abs(value-target),precision),
                'absolute_integral':mp.nstr(absolute,precision),
                'cancellation_ratio':mp.nstr(absolute/abs(value),precision) if value else None,
                'exact':str(exact),'N':N,'p':p,'chords':k,'remaining_sites':L,'internals':internals,
                'cutoff_t':T,'precision':precision}

def main():
    if OUT.exists():raise FileExistsError(OUT)
    start=time.perf_counter();rows=[]
    for p in (16,32):
        for r in range(5):
            ints=tuple(1+i%2 for i in range(r))
            low=calculate(p,ints,60);high=calculate(p,ints,90)
            with mp.workdps(100):
                diff=abs(mp.mpf(low['value'])-mp.mpf(high['value']))
                assert diff<mp.mpf('1e-55')
                assert mp.mpf(high['absolute_error'])<mp.mpf(high['analytic_cutoff_tail'])+mp.mpf('1e-80')
                high['precision_difference']=mp.nstr(diff,95)
            high['lower_precision']=low;rows.append(high)
        print('Fourier finite systems p',p,flush=True)
    paths=['code/compute_central_ordinary_fourier.py','code/central_ordinary_sources.py']
    result=dict(status='passed',created_at_utc=datetime.now(timezone.utc).isoformat(),
        parameters={'p':[16,32],'N':'4p^2','internal_source_count':[0,1,2,3,4],
                    'precisions':[60,90],'cutoff_t':16,'quadrature_breaks':[0,1,2,4,8,16]},
        rows=rows,environment={'python':platform.python_version(),'mpmath':mp.__version__},
        input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        runtime_seconds=time.perf_counter()-start,
        scope='Independent numerical Fourier check against exact rational finite-support coefficients. Analytic tail bounds the omitted t interval only; quadrature and floating roundoff are checked by precision comparison, not interval-certified. No full finite-N many-body simulation.')
    with OUT.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({'status':'passed','rows':len(rows),'seconds':result['runtime_seconds']}),flush=True)

if __name__=='__main__':main()
