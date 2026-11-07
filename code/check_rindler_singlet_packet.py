"""Independent radial quadrature and oscillator checks of the contour formula."""
from pathlib import Path
import json,hashlib
import mpmath as mp
import numpy as np
from scipy.integrate import simpson
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-rindler-packet'


def main():
    path=OUT/'independent-checks.json'
    if path.exists():raise FileExistsError(path)
    study=json.loads((OUT/'study.json').read_text())
    mp.mp.dps=50
    radial=[]
    # Direct high-precision tanh-sinh integral: no float64 nodes or summation.
    for s,e in [(-8,.1),(-8,.2),(0,.2),(1,.1)]:
        e_mp=mp.mpf(str(e));x=mp.exp(s)
        def integrand(n):
            if not n:return mp.mpf(0)
            F=mp.sqrt(8*mp.pi**2/3*n*(1+n*n)/mp.sinh(mp.pi*n))*mp.exp(-e_mp*n*n/2)
            return 2/mp.pi**2*n*mp.sinh(mp.pi*n)*F*mp.besselk(1j*n,x)
        val=mp.re(mp.quad(integrand,list(range(0,41,2))))
        original=next(v for v in study['radial_refinement'] if v['log_x']==s and v['epsilon']==e)
        diff=abs(val-original['low'])
        radial.append(dict(log_x=s,epsilon=e,value_50dps=str(val),reference=original['low'],absolute_difference=float(diff),
            relative_difference=float(diff/abs(val)),passed=bool(diff<mp.mpf('2e-10'))))
        print('independent radial',s,e,flush=True)
    oscillators=[]
    for nu in [.4,1,2]:
        dim=40;h=nu*np.arange(dim);beta=2*np.pi
        rho=np.exp(-beta*h);rho=rho/rho.sum();d=rho**.25
        a=np.diag(np.sqrt(np.arange(1,dim)),1);A=a+a.T
        for eta in [0,.5,2]:
            At=np.exp(1j*eta*(h[:,None]-h[None,:]))*A
            block=(d[:,None]*At)@(d[:,None]*A)
            val=np.trace(block@block)
            c0=1/np.sinh(np.pi*nu)
            ci=np.cos(nu*(eta+1j*np.pi/2))/np.sinh(np.pi*nu)
            wick=c0*c0+2*(ci*ci).real
            oscillators.append(dict(nu=nu,eta=eta,dimension=dim,matrix_value=[val.real,val.imag],wick_value=wick,
                absolute_difference=float(abs(val-wick)),passed=bool(abs(val-wick)<1e-12)))
    grid=np.load(OUT/'radial-grid.npz');norms=[]
    for j,e in enumerate(grid['epsilons']):
        val=simpson(grid['profiles'][:,j]**2,x=grid['log_x'])
        exact=4*np.sqrt(np.pi)/(3*e**1.5)+2*np.sqrt(np.pi)/e**2.5
        norms.append(dict(epsilon=float(e),finite_box_simpson=float(val),exact_full_norm=float(exact),
            relative_difference=float(abs(val/exact-1)),passed=bool(abs(val/exact-1)<1e-4)))
    # Detailed balance and midpoint continuation, evaluated without a thermal
    # partition function for the continuum wedge algebra.
    kms=[]
    for n in [mp.mpf('.2'),mp.mpf('1'),mp.mpf('3')]:
        sc=8*mp.pi/3*n*(1+n*n)/mp.sinh(mp.pi*n)
        plus=mp.exp(mp.pi*n)*sc;minus=mp.exp(-mp.pi*n)*sc
        expected=16*mp.pi/3*n*(1+n*n)
        err=abs((plus-minus)/expected-1)
        kms.append(dict(nu=str(n),commutator_relative_error=str(err),detailed_balance_ratio=str(minus/plus),passed=bool(err<mp.mpf('1e-45'))))
    records=radial+oscillators+norms+kms
    data=dict(status='passed' if all(v['passed'] for v in records) else 'failed',dps=50,
        quadrature='mpmath tanh-sinh, intervals 0,2,...,40; empirical comparison; finite radial box [-8,3] Simpson, no certified spatial tail',
        radial=radial,oscillators=oscillators,finite_box_norms=norms,kms=kms,
        input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in ['code/check_rindler_singlet_packet.py','results/susskind-rindler-packet/study.json','results/susskind-rindler-packet/radial-grid.npz']})
    path.write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(data,indent=2))
    if data['status']!='passed':raise SystemExit(1)


if __name__=='__main__':main()
