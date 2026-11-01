"""Independent four-variable Gaussian quadrature for <f,A B_epsilon f>.

Integrates the continuous Fourier kernel directly; no periodic FFT or
Hermite Hilbert-space truncation. Gauss order changes are not certificates.
"""
from pathlib import Path
import json,hashlib,platform
import numpy as np
import scipy
from scipy.special import roots_hermitenorm
from unitary_flip import multiplier

def calculate(n,epsilon,flat=False,unit=False):
    nodes,weights=roots_hermitenorm(n);weights/=np.sqrt(2*np.pi)
    x=np.sqrt(np.pi)*nodes;p=np.sqrt(np.pi/(1+2*np.pi*epsilon))*nodes
    if unit:a=b=np.ones((n,n),complex);tails=[]
    else:
        a,ta=multiplier(x[:,None],x[None,:],1,flat)
        b,tb=multiplier(p[:,None],p[None,:],1,flat);tails=[ta,tb]
    E=np.exp(1j*x[:,None]*p[None,:]/np.pi)
    inner=E@(b.T*weights[:,None]*weights[None,:])@E.T
    value=2/(1+2*np.pi*epsilon)*np.sum(weights[:,None]*weights[None,:]*a*inner)
    return dict(n=n,epsilon=epsilon,flat=flat,unit=unit,value=[float(value.real),float(value.imag)],
                product_tails=tails,exact_if_unit=1/(1+np.pi*epsilon) if unit else None)

def main():
    root=Path(__file__).resolve().parents[1];out=root/'results/stage4/heat-matrix/direct-left.json'
    if out.exists():raise RuntimeError('Do not overwrite')
    rows=[]
    for n in [128,256,512,768]:
        for e in [.01,1e-10]:
            rows.append(calculate(n,e));print(n,e,rows[-1]['value'],flush=True)
    for e in [.01,1e-10]:
        rows.append(calculate(128,e,unit=True));rows.append(calculate(256,e,flat=True))
    data=dict(rows=rows,versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
        scope='Direct continuous-kernel quadrature of the left word only; empirical Gauss convergence, no rigorous quadrature bound.',
        hashes={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in ['code/direct_heat_left.py','code/unitary_flip.py']})
    out.write_text(json.dumps(data,indent=2))
if __name__=='__main__':main()
