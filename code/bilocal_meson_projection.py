"""Itakura's known bilocal-to-meson transform, applied to saved FLZ modes.

Sources are light-front relative wavefunctions, NOT identified with SYK time
derivatives. Full source moments are analytic; omitted modes stay explicit.
"""
from pathlib import Path
import hashlib,json,platform
from fractions import Fraction
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import eval_chebyu
import scipy

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-bilocal-interface'

def main():
    target=OUT/'projection.json'
    if target.exists():raise FileExistsError(target)
    source=ROOT/'results/susskind-bulk-target/wavefunctions.npz'
    with np.load(source) as saved:
        coeff=saved['coefficients'];mass=np.pi*saved['twice_lambda']
    size=len(coeff);n=np.arange(1,size+1);rows=[]
    # Independent integration grid for the physical source overlaps.
    theta,wt=leggauss(400);theta=(theta+1)*np.pi/2;wt=wt*np.pi/2
    waves=np.sin(theta[:,None]*n)@coeff
    z,wz=leggauss(20)
    for k in range(1,5):
        B=np.zeros(size);same=(k+n)%2==0
        B[same]=1/(1-(k-n[same])**2)-1/(1-(k+n[same])**2)
        bkk=1+1/(4*k*k-1);gkk=bkk/2
        amps=B@coeff/(2*np.sqrt(gkk))
        numerical=(wt*np.sin(theta)*np.sin(k*theta))@waves/(2*np.sqrt(gkk))
        assert np.max(abs(amps-numerical))<2e-12
        probs=abs(amps)**2
        Ckk=float(np.sum(wz*eval_chebyu(k-1,z)**2))
        mean=np.pi*k/bkk
        second=4*k*k*Ckk/bkk
        missing=dict(norm=float(1-probs.sum()),
                     mass2_over_tau=float(mean-probs@mass),
                     mass4_over_tau2=float(second-probs@(mass**2)))
        assert min(missing.values())>-1e-11
        rows.append(dict(source_sine_index=k,first_six_amplitudes=amps.tolist(),first_six_weights=probs.tolist(),
            analytic_mean_mass2_over_tau=float(mean),analytic_second_mass2_over_tau2=float(second),
            mass_moment_ratio=float(second/mean**2),omitted_modes=missing,
            quadrature_max_absolute_difference=float(np.max(abs(amps-numerical)))))
    inputs=['code/bilocal_meson_projection.py','results/susskind-bulk-target/wavefunctions.npz']
    result=dict(status='passed',sources=rows,
        conventions='phi_r normalized in dx; f_k=sin(k theta)/sqrt(B_kk/2), 2x-1=cos(theta); source scale tau=1.',
        environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
        error_scope='Finite 256-mode bulk eigenfunctions with independent quadrature. Only six output modes retained; full source moments and omitted contributions are separately recorded. Not a DSSYK operator map.',
        input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs})
    target.write_text(json.dumps(result,indent=2)+'\n')
    for r in rows:print(r['source_sine_index'],r['first_six_weights'],r['mass_moment_ratio'],r['omitted_modes'])

if __name__=='__main__':main()
