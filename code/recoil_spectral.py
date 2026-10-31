"""Two-point calibration and a declared n=1 recoil-kernel completion.

The finite-mass insertion factor n is not determined by CLP. Model
outputs below are not claims of a unique physical observer correlator.
"""
from pathlib import Path
from datetime import datetime,timezone
import json,hashlib,platform,sys
import mpmath as mp
from recoil_poles import direct_c_pair
ROOT=Path(__file__).resolve().parents[1]

def density(m): return m*m*mp.exp(-2*mp.pi*m)
def spectrum(w,nu): return 1/(8*mp.pi*(mp.cosh(mp.pi*w)+mp.cosh(mp.pi*nu)))
def lr(t,nu):
    return (nu if t==0 else mp.sin(nu*t)/mp.sinh(t))/(4*mp.pi*mp.sinh(mp.pi*nu))
def correction(M,w,nu):
    A=mp.cosh(2*mp.pi*M);B=mp.cosh(mp.pi*w);C=mp.cosh(mp.pi*nu)
    return (1-w*w/(4*M*M))*(A-B)/(A+C)
def insertion_scale(nu): return nu*mp.exp(-mp.pi*nu)/(16*mp.pi**2*mp.sinh(mp.pi*nu))
def joint_kernel(m1,m2,nu,lam=0):
    M=(m1+m2)/2;w=m1-m2
    return density(M)*spectrum(w,nu)*correction(M,w,nu)*(1+lam*w*w/(4*M*M))
def initial_spectrum(m,w,nu,lam=0):return joint_kernel(m,m-w,nu,lam)/density(m)
def partition(m0):
    b=2*mp.pi
    return mp.exp(-b*m0)*(m0*m0/b+2*m0/b**2+2/b**3)
def trace_c(m0):
    b=2*mp.pi
    return (m0+1/b)/(m0*m0+2*m0/b+2/b**2)
def trace_baseline(w,nu,m0):
    w=abs(w)
    return spectrum(w,nu)*mp.exp(-mp.pi*w)*(1+trace_c(m0)*w)
def trace_series(w,nu,m0,terms=20):
    """Exact convergent expansion for n=1, m0>nu/2; truncation recorded."""
    w=abs(w);L=m0+w/2;beta=2*mp.pi
    ss=mp.mpf(0)
    for n in range(terms):
        b=(n+2)*beta
        c=(-1)**n*mp.sinh((n+1)*mp.pi*nu)/mp.sinh(mp.pi*nu)
        pol=(m0*m0+m0*w)/b+(2*m0+w)/b**2+2/b**3
        ss+=c*mp.exp(-b*L)*pol
    return trace_baseline(w,nu,m0)-ss/(4*mp.pi*partition(m0))
def trace_direct(w,nu,m0):
    lo=m0+abs(w)/2
    return spectrum(w,nu)/partition(m0)*mp.quad(lambda M:density(M)*correction(M,w,nu),[lo,lo+1,lo+4,mp.inf])
def trace_tail(nu,m0):return spectrum(0,nu)*density(m0)*correction(m0,0,nu)/partition(m0)
def trace_time_correction(t,nu,m0,terms=20):
    ss=mp.mpf(0);beta=2*mp.pi
    for n in range(terms):
        b=(n+2)*beta;l=b/2
        c=(-1)**n*mp.sinh((n+1)*mp.pi*nu)/mp.sinh(mp.pi*nu)
        A=m0*m0/b+2*m0/b**2+2/b**3;B=m0/b+1/b**2
        integral=2*A*l/(l*l+t*t)+2*B*(l*l-t*t)/(l*l+t*t)**2
        ss+=c*mp.exp(-b*m0)*integral
    return ss/(4*mp.pi*partition(m0))
def trace_series_time_bound(nu,m0,terms=20):
    beta=2*mp.pi;rho=mp.exp(mp.pi*nu-beta*m0);b=(terms+2)*beta
    assert rho<1
    polynomial=2*m0*m0/b**2+8*m0/b**3+8/b**4
    return (mp.exp(-2*beta*m0)*rho**terms*polynomial/
            (2*mp.pi*partition(m0)*(1-mp.exp(-2*mp.pi*nu))*(1-rho)))
def trace_time(t,nu,m0,order=28,dps=45,cutoff=12):
    with mp.workdps(dps):
        t,nu,m0=map(mp.mpf,(t,nu,m0));D=mp.mpf(cutoff);step=mp.mpf('.1')
        nodes,weights=mp.gauss_quadrature(order,'legendre');total=mp.mpf(0)
        segments=int(mp.ceil(D/step))
        for j in range(segments):
            a=D*j/segments;b=D*(j+1)/segments;mid=(a+b)/2;half=(b-a)/2
            total+=half*mp.fsum(weights[k]*2*mp.cos(t*(mid+half*nodes[k]))*
                trace_baseline(mid+half*nodes[k],nu,m0) for k in range(order))
        total-=trace_time_correction(t,nu,m0)
        beta=2*mp.pi;c=trace_c(m0)
        tail=mp.exp(-beta*D)/(2*mp.pi)*((1+c*D)/beta+c/beta**2)
        return total,tail

def run():
    mp.mp.dps=55;checks=[];shapes=[];protocols=[];fourier=[];traces=[];times=[];freedoms=[]
    def check(name,err,tol):
        checks.append(dict(name=name,error=str(abs(err)),tolerance=str(tol),passed=bool(abs(err)<=mp.mpf(tol))))
    nu=mp.mpf('.6')
    for M in map(mp.mpf,[2,5,20]):
        for w in [mp.mpf(0),mp.mpf('.7'),M]:
            m1=M+w/2;m2=M-w/2
            assembled=insertion_scale(nu)*density(m1)*density(m2)*direct_c_pair([m1,m2,nu])/density(M)
            compact=spectrum(w,nu)*correction(M,w,nu)
            check(f'assembled-spectrum-{M}-{w}',assembled/compact-1,'1e-45')
            shapes.append(dict(M=str(M),omega=str(w),assembled=str(assembled),compact=str(compact),
                R=str(correction(M,w,nu)),relative_to_QFT=str(compact/spectrum(w,nu))))
    for m in map(mp.mpf,[2,5,20]):
        for w in map(mp.mpf,['-.7','.7']):
            val=initial_spectrum(m,w,nu)
            reverse=initial_spectrum(m-w,-w,nu)
            check(f'microcanonical-reversal-{m}-{w}',density(m)*val/(density(m-w)*reverse)-1,'1e-45')
            M=m-w/2;A=mp.cosh(2*mp.pi*M);B=mp.cosh(mp.pi*w);C=mp.cosh(mp.pi*nu)
            compact=mp.exp(mp.pi*w)*spectrum(w,nu)*(1-w/m)*(A-B)/(A+C)
            check(f'initial-compact-{m}-{w}',val/compact-1,'1e-45')
            protocols.append(dict(m=str(m),omega=str(w),forward=str(val),reverse_other_state=str(reverse),
                density_initial=str(density(m)),density_final=str(density(m-w)),
                ratio_to_Wightman=str(val/(mp.exp(mp.pi*w)*spectrum(w,nu))),
                same_state_reverse=str(initial_spectrum(m,-w,nu))))
    for t in map(mp.mpf,[0,1,3]):
        val=2*mp.quad(lambda w:mp.cos(w*t)*spectrum(w,nu),[0,1,4,mp.inf])
        check(f'LR-Fourier-{t}',val-lr(t,nu),'1e-40')
        fourier.append(dict(kind='LR',t=str(t),integral=str(val),closed=str(lr(t,nu))))
    # Interior strip: this is an ordinary absolutely convergent integral.
    for eps in map(mp.mpf,['1.2','.7']):
        t=mp.mpf('.8');y=mp.pi-eps
        val=mp.quad(lambda w:mp.exp((-1j*t+y)*w)*spectrum(w,nu),[-mp.inf,-4,0,4,mp.inf])
        closed=lr(t+1j*y,nu)
        check(f'Wightman-interior-strip-{eps}',val-closed,'1e-35')
        fourier.append(dict(kind='regulated-Wightman',epsilon=str(eps),t=str(t),integral=str(val),closed=str(closed)))
    for w in map(mp.mpf,['.2','.7','2']):
        f=lambda x:mp.exp(mp.pi*x)*spectrum(x,nu)
        check(f'QFT-detailed-balance-{w}',f(-w)/f(w)-mp.exp(-2*mp.pi*w),'1e-45')
    for M in map(mp.mpf,[5,20,80]):
        w=mp.mpf('.7')
        for lam in [0,1,2]:
            ratio=joint_kernel(M+w/2,M-w/2,nu,lam)/(density(M)*spectrum(w,nu))
            freedoms.append(dict(M=str(M),omega=str(w),lambda_value=lam,ratio=str(ratio),
                scaled_correction=str(4*M*M*(ratio-1)/(w*w)),predicted_limit=lam-1))
    for m0 in map(mp.mpf,[1,3]):
        for w in map(mp.mpf,[0,1,4]):
            direct=trace_direct(w,nu,m0);series=trace_series(w,nu,m0,24)
            check(f'trace-mass-integral-{m0}-{w}',direct/series-1,'1e-35')
            check(f'trace-series-refinement-{m0}-{w}',trace_series(w,nu,m0,20)/series-1,'1e-32')
            leading=trace_baseline(w,nu,m0)
            bound=2*(1+mp.cosh(mp.pi*nu))*mp.exp(-2*mp.pi*m0)
            check(f'trace-exponential-bound-{m0}-{w}',abs(direct/leading-1),bound)
            traces.append(dict(m0=str(m0),omega=str(w),direct=str(direct),series=str(series),
                leading=str(leading),relative_leading_error=str(abs(direct/leading-1)),bound=str(bound)))
        # One-sided derivative is necessary: the full spectrum has a cusp.
        derivative=mp.diff(lambda w:trace_series(w,nu,m0),mp.mpf(0),direction=1)
        target=-trace_tail(nu,m0)/2
        check(f'trace-cusp-{m0}',derivative/target-1,'1e-30')
        for ti in [8,16,32,64,128]:
            lo,tail=trace_time(ti,nu,m0,28,45);hi,tailhi=trace_time(ti,nu,m0,40,65)
            check(f'trace-time-refinement-{m0}-{ti}',abs(lo/hi-1),'1e-22')
            series_bound=trace_series_time_bound(nu,m0)
            check(f'trace-analytic-truncation-{m0}-{ti}',(tailhi+series_bound)/abs(hi),'1e-24')
            extended=None
            if ti==128:
                extended,_=trace_time(ti,nu,m0,40,65,14)
                check(f'trace-cutoff-extension-{m0}',abs(extended-hi),tailhi*mp.mpf('1.01'))
            times.append(dict(m0=str(m0),t=ti,low=str(lo),value=str(hi),order=[28,40],dps=[45,65],
                relative_refinement=str(abs(lo/hi-1)),baseline_tail_bound=str(tailhi),cutoff=12,step='.1',
                tail_coefficient=str(trace_tail(nu,m0)),tail_ratio=str(hi*ti*ti/trace_tail(nu,m0)),
                series_terms=20,series_time_bound=str(series_bound),extended_value=str(extended) if extended is not None else None,
                relative_analytic_truncation_bound=str((tailhi+series_bound)/abs(hi)),
                scope='n=1 declared completion; baseline Fourier quadrature plus convergent rational correction'))
            print('trace',m0,ti,'done',flush=True)
    inputs=['code/recoil_spectral.py','code/recoil_poles.py','code/recoil_internal.py','code/recoil_vertex.py',
        'sources/chandra-liska-post-2610.30247v1.html','sources/chen-stanford-tang-yang-2607.14042v1.html']
    data=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),status='passed' if all(c['passed'] for c in checks) else 'failed',
        scope='two-point calibration, protocol conversion, normalization freedom and n=1 trace model; no unique finite-mass physical kernel or four-point sewing',
        environment=dict(python=sys.version,platform=platform.platform(),mpmath=mp.__version__,dps=55),nu=str(nu),
        checks=checks,shape_rows=shapes,protocol_rows=protocols,fourier_rows=fourier,freedom_rows=freedoms,
        trace_rows=traces,time_rows=times,hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs})
    out=ROOT/'results/stage4/recoil-spectral';out.mkdir(parents=True,exist_ok=True)
    dest=out/'diagnostics.json'
    assert not dest.exists() or (out/'initial-diagnostics.json').exists()
    dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=data['status'],checks=len(checks),failed=[c for c in checks if not c['passed']]),indent=2))
    if data['status']!='passed':raise SystemExit(1)

if __name__=='__main__':run()
