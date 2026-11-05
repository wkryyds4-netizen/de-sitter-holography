"""One-shot finite-N calibration; exact covariance sums and same-sample singlets."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import platform
import numpy as np
from complex_syk_dictionary import (fermions, creators, covariance, exact_moments,
    sample_hamiltonian, comm, a2_operator, relative, charge_spectrum)

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT/'results/susskind-dictionary'


def main():
    OUT.mkdir(exist_ok=True)
    if (OUT/'calibration.json').exists() or list(OUT.glob('sample-*.npz')):
        raise FileExistsError('Existing calibration retained; choose a new study path.')
    checks = []
    def check(name, error, tol=2e-11):
        checks.append(dict(name=name, error=float(error), tolerance=tol, passed=bool(error<=tol)))
    covariance_rows = []
    for n,p in [(2,1),(3,2),(4,2),(5,2),(5,3)]:
        cs = fermions(n)
        basis = creators(cs,p)
        d = 1<<n
        hsum, asum = 0., 0.
        for a in basis:
            for b in basis:
                x = 2**p * (a@b.conj().T)
                hsum += np.vdot(x,x).real/d
                asum += sum(np.vdot(comm(x,c),comm(x,c)).real/d for c in cs)
        w,a = exact_moments(n,p)
        direct = covariance(n,p)*hsum
        direct_a = covariance(n,p)*asum
        check(f'covariance-H2-{n}-{p}', abs(direct-float(w))/max(1,abs(float(w))))
        check(f'covariance-A2-{n}-{p}', abs(direct_a-float(a))/max(1,abs(float(a))))
        covariance_rows.append(dict(N=n,p=p,H2_formula=float(w),H2_matrix=direct,
                                    A2_formula=float(a),A2_matrix=direct_a))
    samples = []
    for n,p in [(4,1),(3,2),(4,2),(5,2),(6,2),(6,3)]:
        for seed in [1105,1106,1107]:
            cs,j,h = sample_hamiltonian(n,p,seed)
            d = len(h)
            q = sum(c.conj().T@c for c in cs)
            a2 = a2_operator(h,cs)
            o = a2-np.trace(a2)/d*np.eye(d)
            m1 = sum(c.conj().T@(1j*comm(h,c)) for c in cs)
            m2 = -sum(c.conj().T@comm(h,comm(h,c)) for c in cs)
            h2 = h@h
            block_form = sum(c.conj().T@h2@c for c in cs) - (q-2*p*np.eye(d))@h2
            tag = f'{n}-{p}-{seed}'
            check('Hermiticity-'+tag,relative(h,h.conj().T))
            check('charge-'+tag,relative(comm(h,q),np.zeros_like(h)))
            check('M1-'+tag,relative(m1,-1j*p*h))
            check('A2-M2-'+tag,relative(a2,-m2))
            check('charge-block-form-'+tag,relative(a2,block_form))
            inds = [s for s in range(d) if s.bit_count()==p]
            check('lowest-active-sector-'+tag,relative(a2[np.ix_(inds,inds)],p*h2[np.ix_(inds,inds)]))
            om,weight,ch,min_gap = charge_spectrum(h,o,n)
            g = np.vdot(o,o).real/d
            k = np.vdot(comm(h,o),comm(h,o)).real/d
            check('spectral-weight-'+tag,abs(weight.sum()-g)/max(1,g))
            check('spectral-second-moment-'+tag,abs((weight*om**2).sum()-k)/max(1,k))
            cutoff = 1e-10*max(1,np.linalg.norm(h,2))
            static = weight[np.abs(om)<=cutoff].sum()
            dyn = weight[np.abs(om)>cutoff].sum()
            static_small = weight[np.abs(om)<=cutoff/10].sum()
            check('zero-bin-stability-'+tag,abs(static-static_small)/max(1,g))
            if p==1 or n<=p+1:
                check('conserved-control-'+tag,k/max(1,g))
            name = f'sample-{tag}.npz'
            np.savez_compressed(OUT/name,J=j,H=h,A2=a2,omega=om,weight=weight,charge=ch)
            samples.append(dict(N=n,p=p,seed=seed,lambda_value=p*p/n,
                file=name,sha256=hashlib.sha256((OUT/name).read_bytes()).hexdigest(),
                trace_H2=float(np.trace(h2).real/d),G=float(g),K=float(k),
                static_weight=float(static),dynamic_weight=float(dyn),
                dynamic_fraction=float(dyn/g) if g else 0.,zero_cutoff=cutoff,
                min_resolved_gap=min_gap if np.isfinite(min_gap) else None))
    asymptotics = []
    for lam, factor in [(0.25,4),(1.,1),(4.,0.25)]:
        for p in [10,20,40,80]:
            n = int(factor*p*p)
            w,a = exact_moments(n,p)
            asymptotics.append(dict(N=n,p=p,lambda_value=lam,
                                   lambda_H2=float(w)*lam,
                                   p_fundamental_mu2=2*p*float(a)/n))
    report = dict(created_at_utc=datetime.now(timezone.utc).isoformat(),
        state='normalized full trace I/2^N, beta=mu=0', Js=1.,
        covariance='Hermitian complex Gaussian: E[J_IJ J_KL]=v delta_IL delta_JK; diagonal real variance v, offdiagonal E|J|^2=v',
        time='string convention with H_s; H_c=p H_s is a separate time-unit conversion',
        python=platform.python_version(),numpy=np.__version__,
        scope='Finite-N calibration and sample spectra only; no asymptotic singlet spectrum or bulk-mass claim.',
        code_sha256={str(x.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(x.read_bytes()).hexdigest()
                     for x in [Path(__file__),ROOT/'code/complex_syk_dictionary.py']},
        covariance_rows=covariance_rows,samples=samples,asymptotics=asymptotics,
        checks=checks,passed=all(c['passed'] for c in checks))
    with (OUT/'calibration.json').open('x',encoding='utf-8') as f:
        json.dump(report,f,indent=2);f.write('\n')
    print(json.dumps({'passed':report['passed'],'checks':len(checks),'samples':len(samples),
                      'failed':[c for c in checks if not c['passed']]}))
    if not report['passed']:
        raise SystemExit(1)


if __name__=='__main__':
    main()
