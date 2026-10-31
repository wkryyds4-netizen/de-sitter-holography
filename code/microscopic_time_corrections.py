"""C.39 time corrections and the complete conformal-antipodal C.6 kernel.

All statements are within Narovlansky 2506.02109v2's semiclassical integral.
No finite-observer mass dictionary, finite-G match, or finite-lambda completion.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import platform
import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
SOURCE = 'https://arxiv.org/html/2506.02109v2'
HALF = mp.mpf('0.5')


def barnes_i(z):
    """Integral of Gamma(a+z+iw)Gamma(b+z+iw)Gamma(a-iw)Gamma(b-iw)."""
    return 2*mp.pi*mp.gamma(z+HALF)*mp.gamma(z+1)**2*mp.gamma(z+mp.mpf('1.5'))/mp.gamma(2*z+2)


def barnes_d(n):
    return mp.digamma(n+HALF)+mp.digamma(n+1)-2*mp.log(2)


def gamma_weight(w, n):
    a, b = mp.mpf('.25'), mp.mpf('.75')
    return mp.gamma(a+n+1j*w)*mp.gamma(b+n+1j*w)*mp.gamma(a-1j*w)*mp.gamma(b-1j*w)


def direct_i_k(n):
    intervals = [0, 1, 3, 8, 16, mp.inf]
    iv = 2*mp.re(mp.quad(lambda w: gamma_weight(w,n), intervals))
    kv = 2*mp.re(mp.quad(lambda w: gamma_weight(w,n)*(mp.digamma(n+mp.mpf('.25')+1j*w)+mp.digamma(n+mp.mpf('.75')+1j*w)), intervals))
    return iv, kv


def b_n(n, s0, t):
    return mp.log(8*s0)-t-mp.digamma(n+HALF)


def c39_term(n, s0, t):
    return mp.rf(HALF,n)**2*(mp.exp(2*t)/(64*s0*s0))**n*b_n(n,s0,t)/mp.log(s0)


def mixed_barnes(w):
    return 2*mp.pi*mp.gamma(HF-w*1j)*mp.gamma(mp.mpf('1.5')+w*1j)


HF = HALF


def raw_combined_kernel(w, s0):
    """Sum all four C.6 sectors before taking w->0; excludes exp(-2itw)."""
    big_l = mp.log(2*s0)
    first = mp.exp(2j*big_l*w)*mp.gamma(1j*w)**2*barnes_i(-1j*w)**2
    fourth = mp.exp(-2j*big_l*w)*mp.gamma(-1j*w)**2*barnes_i(1j*w)**2
    mix = mp.gamma(1j*w)*mp.gamma(-1j*w)*(mixed_barnes(w)**2+mixed_barnes(-w)**2)
    return first+fourth+mix


def reduced_kernel(w, s0):
    """raw_combined_kernel/(2*pi^6), evaluated with removable cancellation."""
    ell = mp.log(8*s0)
    if not w:
        return (2*(ell-mp.digamma(HALF))**2-4-mp.pi**2/3)/mp.pi**2
    if abs(w)<mp.mpf('1e-8'):
        w2=w*w
        a_minus_one=(-(4+mp.pi**2/3)*w2+(4*mp.pi**2/3+2*mp.pi**4/15)*w2**2
                     +(-8*mp.pi**4/15-17*mp.pi**6/315)*w2**3
                     +(68*mp.pi**6/315+62*mp.pi**8/2835)*w2**4)
    else:
        a_minus_one=(1-4*w*w)*mp.tanh(mp.pi*w)/(mp.pi*w)-1
    half_phase=w*ell-mp.im(mp.loggamma(HALF+1j*w))
    return (a_minus_one+2*mp.sin(half_phase)**2)/(mp.sinh(mp.pi*w)**2*mp.cosh(mp.pi*w))


def full_c6_ratio(s0, t, refined=False):
    """R_micro/C0, with the same fixed projection scale in both quantities."""
    points = [0, mp.mpf('.125'),mp.mpf('.25'),mp.mpf('.5'),1,2,4,8,16,mp.inf]
    if refined:
        points=[mp.mpf(k)/8 for k in range(17)]+[3,4,6,8,12,16,mp.inf]
    return mp.pi/mp.log(s0)*mp.quad(lambda w: mp.cos(2*t*w)*reduced_kernel(w,s0), points)


def main():
    mp.mp.dps=50
    checks=[]
    def num(x): return mp.nstr(x,40)
    def record(name, error, tolerance='1e-35'):
        checks.append({'name':name,'error':num(error),'tolerance':tolerance,'passed':bool(error<mp.mpf(tolerance))})
    direct=[]
    for n in [0,1,2,4]:
        iv,kv=direct_i_k(n)
        exact=barnes_i(n)
        record('direct Barnes integral n='+str(n),abs(iv/exact-1))
        record('direct digamma-weighted integral n='+str(n),abs(kv/exact-barnes_d(n)))
        record('derivative of Barnes expression n='+str(n),abs(mp.diff(lambda z:mp.log(barnes_i(z)),mp.mpf(n))-barnes_d(n)))
        direct.append({'n':n,'I_integral':num(iv),'I_closed':num(exact),'K_over_I':num(kv/iv),'D_closed':num(barnes_d(n))})
        print('Barnes and digamma verified n=',n,flush=True)
    s0=mp.exp(50);t=mp.mpf(6);l=mp.log(s0)
    coeffs=[]
    for n in range(9):
        direct_bracket=-2j*mp.digamma(n+1)-mp.pi+2j*(t-1j*mp.pi/2)-2j*mp.log(2*s0)+2j*barnes_d(n)
        raw=(2*mp.pi*1j*mp.exp(-1j*mp.pi*n)*(-1)**n/mp.factorial(n)**2/4**n
             *barnes_i(n)**2*direct_bracket/(4*mp.pi*l*barnes_i(0)**2))
        closed=mp.rf(HALF,n)**2/64**n*b_n(n,s0,t)/l
        record('C.39 bracket and phases n='+str(n),abs(raw-closed))
        harmonic=2*sum(mp.mpf(1)/(2*k+1) for k in range(n))
        record('odd-harmonic shift n='+str(n),abs(b_n(0,s0,t)-b_n(n,s0,t)-harmonic))
        coeffs.append({'n':n,'coefficient_after_factoring_exp2nt_over_s02n':num(closed),
                       'B_n':num(b_n(n,s0,t)),'shift_from_B0':num(harmonic)})
    for wstr in ['0.01','0.1','0.7','2']:
        w=mp.mpf(wstr)
        raw=raw_combined_kernel(w,s0)
        alt=2*mp.pi**6*reduced_kernel(w,s0)
        record('four-sector kernel reduction w='+wstr,abs(raw-alt)/max(1,abs(alt)))
    # Independently integrate the mixed Barnes factor, whose shifts differ from I.
    for wstr in ['0.1','0.7']:
        w=mp.mpf(wstr); a=mp.mpf('.25');b=mp.mpf('.75')
        mixed=mp.quad(lambda x:mp.gamma(a-1j*w+1j*x)*mp.gamma(b+1j*x)*mp.gamma(a-1j*x)*mp.gamma(b+1j*w-1j*x),[-mp.inf,-8,-3,-1,0,1,3,8,mp.inf])
        record('direct mixed Barnes factor w='+wstr,abs(mixed/mixed_barnes(w)-1))
    # An explicit approach to the same antipodal point from the negative strip.
    eps_rows=[];previous=None
    for epsstr in ['0.02','0.002','0.0002']:
        eps=mp.mpf(epsstr);n=1
        iv=mp.quad(lambda w:mp.exp(2*eps*w)*gamma_weight(w,n),[-mp.inf,-8,-3,-1,0,1,3,8,mp.inf])
        kv=mp.quad(lambda w:mp.exp(2*eps*w)*gamma_weight(w,n)*(mp.digamma(n+mp.mpf('.25')+1j*w)+mp.digamma(n+mp.mpf('.75')+1j*w)),[-mp.inf,-8,-3,-1,0,1,3,8,mp.inf])
        bracket=-eps+2j*(t-mp.log(2*s0)-mp.digamma(n+1))+2j*mp.re(kv/iv)
        value=2*mp.pi*1j*mp.exp(1j*n*eps)*abs(iv)**2*bracket/(mp.factorial(n)**2*4**n*4*mp.pi*l*barnes_i(0)**2)
        target=mp.rf(HALF,n)**2/64**n*b_n(n,s0,t)/l
        error=abs(value-target)
        if previous is not None: record('negative-strip convergence epsilon='+epsstr,abs(error/previous-mp.mpf('.1')),'1e-3')
        eps_rows.append({'epsilon':epsstr,'coefficient_real':num(mp.re(value)),'coefficient_imag':num(mp.im(value)),'error_to_boundary':num(error)})
        previous=error
    print('all coefficient, kernel and regulator checks complete',flush=True)
    full=[]
    for logs,tstr in [(20,'5'),(20,'10'),(50,'6')]:
        sv=mp.exp(logs);tv=mp.mpf(tstr)
        with mp.workdps(40):v40=full_c6_ratio(sv,tv)
        v50=full_c6_ratio(sv,tv,True)
        record('full kernel precision and subdivision log_s='+str(logs)+',t='+tstr,abs(v50-v40),'1e-28')
        n0=c39_term(0,sv,tv);n1=c39_term(1,sv,tv);n2=c39_term(2,sv,tv)
        full.append({'log_s0':logs,'t':tstr,'full_C6_over_C0':num(v50),'C39_n0':num(n0),'C39_n1':num(n1),'C39_n2':num(n2),
                     'full_minus_C39_through_n2':num(v50-n0-n1-n2),'precision_difference':num(abs(v50-v40))})
        print('full C.6 log_s=',logs,'t=',tstr,'ratio=',num(v50),flush=True)
    failed=sum(not row['passed'] for row in checks)
    files=[Path(__file__),ROOT/'sources/narovlansky-2506.02109v2.html']
    result={'generated_utc':datetime.now(timezone.utc).isoformat(),'status':'failed' if failed else 'passed','checks':len(checks),'failed':failed,
            'check_details':checks,'precision_dps':50,'source':SOURCE,'source_equations':['C.6','C.9','C.13','C.18','C.21','C.39','C.40'],
            'normalization':'R_micro/C0; C0=mathfrak_J log(s0)/(4*pi^2*T), fixed same mathfrak_J; t=2*pi*T*t_phys',
            'scope':'Exact reduction of the already-semiclassical C.6 integral at two conformal antipodal probes; C.39 coefficient corrections remain an early-time asymptotic block. No finite-lambda, finite-G or observer-mass match.',
            'coefficient_test_parameters':{'log_s0':'50','t':'6'},'direct_barnes':direct,'C39_coefficients':coeffs,'negative_strip_limit':eps_rows,'full_C6_examples':full,
            'code_sha256':{Path(__file__).name:hashlib.sha256(files[0].read_bytes()).hexdigest()},
            'input_sha256':{files[1].relative_to(ROOT).as_posix():hashlib.sha256(files[1].read_bytes()).hexdigest()},
            'environment':{'python':platform.python_version(),'mpmath':mp.__version__}}
    out=ROOT/'results/stage4/time-corrections.json'
    out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(f'{len(checks)-failed} passed, {failed} failed',flush=True)
    if failed:
        for c in checks:
            if not c['passed']:print(c)
        raise SystemExit(1)


if __name__=='__main__':main()
