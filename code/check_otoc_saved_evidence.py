"""Replay saved arithmetic, high-precision cutoff bounds, and phase diagnostics."""
from pathlib import Path
import hashlib,json,platform
import numpy as np
import mpmath as mp
from singlet_otoc import Kernel

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-singlet-otoc'


def main():
    target=OUT/'error-replay.json'
    if target.exists():raise FileExistsError(target)
    data=json.loads((OUT/'otoc-numerics.json').read_text());mp.mp.dps=70
    bounds=[];phases=[];count=0
    for g in data['groups']+data['scrambling_scaling']:
        lam=mp.mpf(str(g['lambda_value']));q=mp.exp(-4*lam);r=q*q;chi=q**4
        c=g['cutoff'];L=c['L'];M=c['inner'];tmax=max(v['time'] for v in g['records'])
        a=2*mp.mpf(str(tmax))/mp.sqrt(lam*(1-q));K=1/((1-q)*(1-q**2)*(1-q**3)*(1-r))
        def tail(d):return a**d/mp.factorial(d)/(1-a/(d+1))
        eps0=2*tail(L);delta=eps0+2*tail(L-M)+2*r**M
        bound=(4+2*chi*K)*delta
        assert bound<mp.mpf('1e-10')
        assert abs(float(bound)/c['total_otoc_truncation_bound']-1)<1e-10
        bounds.append(dict(lambda_value=str(lam),L=L,inner=M,tmax=tmax,otoc_bound_70dps=str(bound),
            commutator_bound_70dps=str(2*(2*(1+chi)*eps0+bound))))
        for v in g['records']:
            assert abs(v['otoc']-(2*v['uncrossed']+v['crossed_weighted']))<1e-14
            assert abs(v['crossed_weighted']-float(chi)*v['crossed_reduced'])<1e-14
            assert abs(v['commutator_squared']-2*(v['f22']-v['otoc']))<1e-14
            count+=1
    for g in data['groups']:
        k=Kernel(g['lambda_value'],g['cutoff']['L'])
        for row in g['records']:
            t=row['time'];minus=(k.u*np.exp(-1j*t*k.e)[None,:])@k.u.T
            v=minus.conj()@(k.rn*minus[:,0]);terms=v*v;value=terms.sum()
            assert abs(value.real-row['uncrossed'])<1e-12
            phases.append(dict(lambda_value=g['lambda_value'],time=t,
                uncrossed_phase_condition=float(sum(abs(terms))/max(abs(value),1e-300)),
                total_channel_condition=(2*abs(row['uncrossed'])+abs(row['crossed_weighted']))/max(abs(row['otoc']),1e-300),
                uncrossed_imag=float(value.imag)))
    inputs=['code/check_otoc_saved_evidence.py','code/singlet_otoc.py','results/susskind-singlet-otoc/otoc-numerics.json']
    payload=dict(status='passed',arithmetic_records=count,bound_precision=70,bounds=bounds,phase_records=phases,
        environment=dict(python=platform.python_version(),numpy=np.__version__,mpmath=mp.__version__),
        input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs})
    target.write_text(json.dumps(payload,indent=2)+'\n')
    print(json.dumps(dict(status='passed',arithmetic_records=count,bound_records=len(bounds),phase_records=len(phases),
        max_uncrossed_phase_condition=max(v['uncrossed_phase_condition'] for v in phases),
        max_total_channel_condition=max(v['total_channel_condition'] for v in phases))))

if __name__=='__main__':main()
