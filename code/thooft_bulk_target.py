"""Reproduce the FLZ alpha=0 meson target, independently in two representations.

FLZ 0905.2280 Eqs (1.1)-(1.5), (8.5)-(8.7), and section 5.
No DSSYK-to-spatial-operator map is assumed or fitted.
"""
from pathlib import Path
import hashlib, json, platform
import numpy as np
import scipy
from scipy.linalg import eigh
from numpy.polynomial.legendre import leggauss

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / 'results/susskind-bulk-target'
# First six published 2*lambda values, FLZ Tables 1 and 2.
REFERENCE = np.array([.73706174629269, 1.7537313369175, 2.7481609123706,
                      3.7510575817054, 4.7492953810375, 5.7504926236487])

def chebyshev(size):
    n = np.arange(1, size+1, dtype=float)
    m = n[:, None]
    same = ((m+n) % 2 == 0)
    a, b = 1-(m-n)**2, 1-(m+n)**2
    B = np.zeros((size, size))
    B[same] = (1/a[same] - 1/b[same])
    # B_mn = integral sin(m theta) sin(n theta) sin(theta) dtheta.
    # D c = 2 lambda B c; eigenvalues of K are 1/lambda.
    K = 2*B/np.sqrt(m*n)
    inv, vec = eigh(K, subset_by_index=[size-6, size-1])
    inv, vec = inv[::-1], vec[:, ::-1]
    lam = 1/inv
    residual = np.linalg.norm(K@vec-vec*inv, axis=0)
    # Reconstruct phi=sum c_n sin(n theta), normalize in dx.
    coeff = vec/np.sqrt(n[:, None])
    coeff /= np.sqrt(np.sum(coeff*(B@coeff), axis=0)/2)
    gram = coeff.T@B@coeff/2
    parity = np.sum(coeff*(B@(coeff*((-1)**(n+1))[:,None])), axis=0)/2
    return dict(size=size, twice_lambda=(2*lam).tolist(), mass_squared_over_tension=(2*np.pi*lam).tolist(),
        relative_reference_error=(2*lam/REFERENCE-1).tolist(),
        matrix_residual_max=float(residual.max()), orthogonality_error=float(np.max(abs(gram-np.eye(6)))),
        reflection_parity=parity.tolist()), coeff

def fredholm(limit, size):
    nodes, weights = leggauss(size)
    t, w = limit*nodes, limit*weights
    z = t[:,None]-t[None,:]
    S = np.ones_like(z)
    np.divide(z, np.sinh(z), out=S, where=z!=0)
    f = t/np.tanh(t)
    factors = np.sqrt(w/f)
    K = factors[:,None]*S*factors[None,:]
    inv = eigh(K, subset_by_index=[size-6,size-1], eigvals_only=True)[::-1]
    vals = 2/inv
    return dict(limit=limit, nodes=size, twice_lambda=vals.tolist(),
                relative_reference_error=(vals/REFERENCE-1).tolist())

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    target = OUT/'spectrum.json'
    if target.exists():
        raise FileExistsError(target)
    rows=[]
    for size in [32,64,128,256]:
        row, coeff=chebyshev(size)
        rows.append(row)
        print('Chebyshev',size,'max relative error',max(abs(np.array(row['relative_reference_error']))),flush=True)
    other=[]
    for limit,size in [(16,160),(24,240),(24,360),(32,480),(40,600),(40,800)]:
        row=fredholm(limit,size);other.append(row)
        print('Fredholm',limit,size,'max relative error',max(abs(np.array(row['relative_reference_error']))),flush=True)
    (OUT/'refinement-checkpoint.json').write_text(json.dumps(dict(chebyshev=rows,fredholm=other),indent=2)+'\n')
    assert max(abs(np.array(rows[-1]['relative_reference_error']))) < 2e-11
    assert max(abs(np.array(other[-1]['relative_reference_error']))) < 2e-9
    assert np.max(abs(np.array(rows[-1]['twice_lambda'])/other[-1]['twice_lambda']-1)) < 2e-9
    assert rows[-1]['orthogonality_error'] < 1e-11
    assert np.max(abs(np.array(rows[-1]['reflection_parity'])-(-1.)**np.arange(6)))<1e-11
    grid=np.linspace(.001,.999,999);theta=np.arccos(2*grid-1)
    wave=np.sin(theta[:,None]*np.arange(1,257))@coeff
    for j in range(6):
        if wave[np.argmax(abs(wave[:,j])),j]<0:
            coeff[:,j]*=-1;wave[:,j]*=-1
    np.savez(OUT/'wavefunctions.npz',x=grid,phi=wave,coefficients=coeff,
             twice_lambda=np.array(rows[-1]['twice_lambda']))
    result=dict(status='passed',source='https://arxiv.org/html/0905.2280',
        conventions=dict(alpha=0, FLZ_g_squared='tau=g_YM^2*N', FLZ_lambda='nu, NOT the DSSYK coupling lambda',
                         mass_squared='2*pi*tau*nu', microscopic_limit='candidate planar bulk target only'),
        reference_twice_lambda=REFERENCE.tolist(),chebyshev=rows,fredholm=other,
        error_scope='Matrix residuals, basis/domain refinements and independent published values. No certified continuum eigenvalue enclosure; no microscopic holographic match.',
        environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
        input_sha256={'code/thooft_bulk_target.py':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    target.write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':
    main()
