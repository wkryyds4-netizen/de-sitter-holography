"""Plot candidate boost and translated overlaps; no fit to new SYK data."""
from pathlib import Path
import json,hashlib,platform
import numpy as np
import scipy
from scipy.special import kv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-translations'

def main():
    target=OUT/'plot-data.json'
    if target.exists():raise FileExistsError(target)
    eta=np.linspace(0,3,121);a=np.linspace(0,2,121);curves=[]
    fig,ax=plt.subplots(1,2,figsize=(9,3.3))
    ax[0].plot(eta,1/np.cosh(eta/2)**4,'k',label='all positive masses')
    for m in [.5,1,2]:
        alpha=2+1j*m*a[1:]/2;beta=1j*m*a[1:]/2
        vals=np.r_[1+0j,(16/3)*(beta/alpha)**2*kv(4,2*np.sqrt(alpha*beta))]
        curves.append(dict(m=m,real=vals.real.tolist(),imag=vals.imag.tolist()))
        ax[1].plot(a,vals.real,label=fr'$m={m}$')
    ax[0].set(xlabel=r'boost parameter $\eta$',ylabel=r'$\langle\psi|e^{-i\eta B}|\psi\rangle$')
    ax[1].set(xlabel=r'inertial time displacement $a$',ylabel=r'$\mathrm{Re}\langle\psi|e^{-iaP^0}|\psi\rangle$')
    for x in ax:x.legend(fontsize=8);x.grid(alpha=.2)
    fig.tight_layout();fig.savefig(OUT/'candidate-translations.png',dpi=180);plt.close(fig)
    target.write_text(json.dumps(dict(scope='Candidate one-particle kinematics, d=2, mu=m; not measured microscopic translations',eta=eta.tolist(),a=a.tolist(),curves=curves,
        environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,matplotlib=matplotlib.__version__),
        input_sha256={'code/plot_singlet_translations.py':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}),indent=2)+'\n')
    print('Candidate comparison figure saved.')
if __name__=='__main__':main()
