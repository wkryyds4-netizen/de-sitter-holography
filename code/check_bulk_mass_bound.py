"""Numerical controls for the analytic planar form bound, not SYK error bounds."""
from pathlib import Path
import json, hashlib, platform
import numpy as np
import scipy
from scipy.linalg import eigh
from scipy.special import eval_chebyu
from numpy.polynomial.legendre import leggauss

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-bulk-target'

def main():
    target=OUT/'mass-bound.json'
    if target.exists():raise FileExistsError(target)
    records=[]
    for size in [32,64]:
        n=np.arange(1,size+1);m=n[:,None];same=(m+n)%2==0
        B=np.zeros((size,size));aa=1-(m-n)**2;bb=1-(m+n)**2
        B[same]=1/aa[same]-1/bb[same]
        # C_mn = integral_{-1}^1 U_(m-1)(z) U_(n-1)(z) dz.
        # The Gauss rule is exact for these polynomials in exact arithmetic.
        z,w=leggauss(size+4)
        values=eval_chebyu(n[:,None]-1,z[None,:]);C=(values*w)@values.T
        Q0=np.diag(np.pi**2*n/2);V=2*C;G=B/2
        lo=float(eigh(Q0-V,Q0,eigvals_only=True)[0]);assert lo>-1e-11
        base=eigh(Q0,G,subset_by_index=[0,5],eigvals_only=True)
        for p in [4,8,16,32]:
            alpha=np.pi/p**2
            shifted=eigh(Q0+alpha*V,G,subset_by_index=[0,5],eigvals_only=True)
            relative=shifted/base-1
            assert relative.min()>-1e-11 and relative.max()<=alpha+1e-11
            records.append(dict(size=size,p=p,alpha=alpha,relative_shifts=relative.tolist(),
                form_bound_slack_min=lo,finite_basis_only=True))
    target.write_text(json.dumps(dict(status='passed',records=records,
        scope='Analytic proof is q0<=qa<=(1+a)q0. These are finite-dimensional arithmetic controls and do not estimate eigenvalue errors at nonzero a.',
        environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
        input_sha256={'code/check_bulk_mass_bound.py':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}),indent=2)+'\n')
    print('8 finite-basis mass comparisons passed; full continuum bound follows from the positive form.')

if __name__=='__main__':main()
