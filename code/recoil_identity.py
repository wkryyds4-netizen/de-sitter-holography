"""Identity-limit diagnostics, not a physical dS sewing construction.

Radius one; first fixed-mass kappa -> infinity. Independent real-axis
quadratures test a normalized continued two-point spectrum and pole residues.
"""
from pathlib import Path
import hashlib
import json
import platform
from datetime import datetime, timezone
import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/stage4/recoil-identity'


def density(w, e):
    # Stable near the pinch, including negative e.
    return mp.sin(mp.pi*e)/(4*(1-e)*(mp.sinh(mp.pi*w/2)**2
                                      + mp.sin(mp.pi*e/2)**2))


def characteristic(t, e):
    return mp.mpf(1) if not t else mp.sinh((1-e)*t)/((1-e)*mp.sinh(t))


def recoil(w, e, M):
    A = mp.cosh(2*mp.pi*M)
    return (1-w*w/(4*M*M))*(A-mp.cosh(mp.pi*w))/(A-mp.cos(mp.pi*e))


def integral(e, t=0, power=0, M=None, lower=None, lam=0):
    D = mp.inf if M is None else 2*(M-lower)
    points = sorted(set([mp.mpf(0)] + [x for x in
        [abs(e), 10*abs(e), mp.mpf(1), mp.mpf(4), mp.mpf(12)] if 0 < x < D] + [D]))
    def f(w):
        a = density(w,e)*mp.cos(w*t)*w**power
        if M is not None:
            a *= recoil(w,e,M)*(1+lam*w*w/(4*M*M))
        return a
    return 2*mp.quad(f, points)


def run():
    mp.mp.dps = 50
    rows, checks = [], []
    def check(name, observed, expected, tol):
        error = abs(observed-expected)
        checks.append(dict(name=name, observed=mp.nstr(observed,48),
            expected=mp.nstr(expected,48), absolute_error=mp.nstr(error,12),
            tolerance=str(tol), passed=bool(error <= mp.mpf(tol))))
    for es in ['.2','.05','.01','.002']:
        e = mp.mpf(es)
        for t in [mp.mpf(0),mp.mpf('.7'),mp.mpf(3)]:
            val = integral(e,t)
            check(f'positive_fourier_e{es}_t{t}',val,characteristic(t,e),'1e-30')
        moment = integral(e,power=2)
        check(f'second_moment_e{es}',moment,(2*e-e*e)/3,'1e-30')
        rows.append(dict(kind='approximate_identity',epsilon=es,
            variance=mp.nstr(moment,40),t3=mp.nstr(characteristic(3,e),40)))
    for ds in ['.2','.05','.01']:
        d=mp.mpf(ds)
        for t in [mp.mpf(0),mp.mpf('.7'),mp.mpf(3)]:
            val=integral(-d,t)
            naive=-(1-d)/(1+d)*characteristic(t,d)
            correction=2*mp.cosh(d*t)/(1+d)
            check(f'negative_real_axis_d{ds}_t{t}',val,naive,'1e-30')
            check(f'crossing_residues_d{ds}_t{t}',val+correction,
                  characteristic(t,-d),'1e-30')
            rows.append(dict(kind='pole_crossing',delta=ds,t=str(t),
                real_axis=mp.nstr(val,40),residue_correction=mp.nstr(correction,40),
                analytically_continued=mp.nstr(characteristic(t,-d),40)))
    finite=[]
    for Ms in ['1','3']:
        M,lower=mp.mpf(Ms),mp.mpf('.2')
        D=2*(M-lower)
        for es in ['.1','.02','.004']:
            e=mp.mpf(es); v=(2*e-e*e)/3
            A=mp.cosh(2*mp.pi*M)
            uniform_bound=v/D**2+v/(4*M*M)+mp.sin(mp.pi*e)/(
                (1-e)*(A-mp.cos(mp.pi*e)))*(D-D**3/(12*M*M))
            val=integral(e,M=M,lower=lower)
            tail=1-integral(e,M=None)+2*mp.quad(lambda w:density(w,e),[D,mp.inf])
            truncated_moment=2*mp.quad(lambda w:w*w*density(w,e),[0,e,10*e,D])
            exact_loss=tail+truncated_moment/(4*M*M)+mp.sin(mp.pi*e)/(
                (1-e)*(A-mp.cos(mp.pi*e)))*(D-D**3/(12*M*M))
            check(f'finite_mass_loss_M{Ms}_e{es}',1-val,exact_loss,'1e-30')
            check(f'finite_mass_bound_M{Ms}_e{es}',max(0,1-val-uniform_bound),0,'1e-30')
            for lam in [0,2]:
                t=mp.mpf(2)
                q=integral(e,t,M=M,lower=lower,lam=lam)
                bound=uniform_bound+lam*v/(4*M*M)
                check(f'uniform_fourier_M{Ms}_e{es}_lambda{lam}',
                      max(0,abs(q-characteristic(t,e))-bound),0,'1e-30')
                with mp.workdps(65):
                    refined=integral(mp.mpf(es),t,M=M,lower=lower,lam=lam)
                check(f'quadrature_refinement_M{Ms}_e{es}_lambda{lam}',q,refined,'1e-38')
                finite.append(dict(M=Ms,lower='.2',epsilon=es,lam=lam,
                    t='2',integral=mp.nstr(q,48),refined=mp.nstr(refined,48),
                    distance_to_identity=mp.nstr(abs(q-1),30),
                    uniform_fourier_error_bound=mp.nstr(bound,30)))
    weights=[]
    for ks in ['2','10','100']:
        k=mp.mpf(ks)
        for nu in [mp.mpf('.6'),mp.j,-mp.j]:
            hs=[]
            for sign in [1,-1]:
                b=mp.sqrt(sign*mp.j/k); Q=b+1/b
                alpha=2*mp.pi*(1-nu/k)
                P=sign*mp.j*alpha/(4*mp.pi*b)
                h=Q*Q/4+P*P
                formula=(1-sign*mp.j*nu)/2+sign*mp.j*(1+nu*nu)/(4*k)
                check(f'weight_k{ks}_nu{nu}_sign{sign}',h,formula,'1e-45')
                vacuum_alpha=2*mp.pi*(1+b*b)
                vacuum_P=sign*mp.j*vacuum_alpha/(4*mp.pi*b)
                check(f'vacuum_k{ks}_nu{nu}_sign{sign}',Q*Q/4+vacuum_P**2,0,'1e-45')
                hs.append(mp.nstr(h,35))
            weights.append(dict(kappa=ks,nu=str(nu),h_plus=hs[0],h_minus=hs[1]))
    # No quadrature claim for the high-time scaling samples.
    scaling=[]
    for es in ['.05','.01','.002']:
        e=mp.mpf(es)
        for tau in [mp.mpf('.5'),mp.mpf(1),mp.mpf(2)]:
            scaling.append(dict(epsilon=es,tau=str(tau),
                exact=mp.nstr(characteristic(tau/e,e),40),limit=mp.nstr(mp.exp(-tau),40)))
    inputs=['code/recoil_identity.py','sources/collier-eberhardt-zhang-2304.13650v3.html',
            'sources/chandra-liska-post-2610.30247v1.html','sources/collier-et-al-2501.01486v1.html',
            'notes/stage4-recoil-spectral.md']
    data=dict(status='passed' if all(x['passed'] for x in checks) else 'failed',
        created_at_utc=datetime.now(timezone.utc).isoformat(),python=platform.python_version(),
        mpmath=mp.__version__,precision=50,refinement_precision=65,checks=checks,
        rows=rows,finite_mass=finite,weights=weights,scaling=scaling,
        scope='Reduced two-point continuation and exact label algebra; not finite-kappa sewing, not interval quadrature certification.',
        hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs})
    OUT.mkdir(parents=True,exist_ok=True)
    target=OUT/'diagnostics.json'
    if target.exists():raise FileExistsError('Keep prior diagnostic run; use a new filename.')
    target.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=data['status'],checks=len(checks),
        failed=[c['name'] for c in checks if not c['passed']],
        max_refinement=max(float(c['absolute_error']) for c in checks if 'refinement' in c['name']))))


if __name__=='__main__':run()
