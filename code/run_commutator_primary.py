"""Finite-lambda mixed commutator and independent unequal-probe transfer checks."""
from pathlib import Path
import json,hashlib,platform,time
import numpy as np
import mpmath as mp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mixed_cubic_otoc import Kernel,cutoff
ROOT=Path(__file__).resolve().parent.parent;OUT=ROOT/'results/susskind-commutator-primary'

def nonsymmetric(lam,t,L):
    q=mp.exp(-4*lam);a=q**2;b=q**3;T=mp.matrix(L)
    for n in range(L-1):T[n+1,n]=1/mp.sqrt(lam);T[n,n+1]=(1-q**(n+1))/((1-q)*mp.sqrt(lam))
    minus=mp.expm(-1j*t*T);plus=minus.conjugate();M=minus*mp.diag([a**j for j in range(L)])*plus
    f=[mp.mpf(1)];g=[mp.mpf(1)]
    for j in range(1,L):f.append(f[-1]*(1-q**j));g.append(g[-1]*(1-b*b*q**(j-1)))
    terms=[]
    for i in range(L):
      for l in range(L-i):
       P=f[l+i]*g[i]/(f[l]*f[i])
       for m in range(L-i):terms.append(P*a**i*b**(l+m)*plus[0,m+i]*M[m,l]*minus[l+i,0])
    reduced=mp.fsum(terms)
    return q**6*reduced,mp.fsum(abs(v) for v in terms)/abs(reduced)
def main():
    target=OUT/'numerics.json'
    if target.exists():raise FileExistsError(target)
    start=time.monotonic();rows=[]
    for lam in [.02,.01,.005]:
        cut=cutoff(lam,.35,2,3);k=Kernel(lam,cut['L']);points=[]
        for t in [0.,.1,.2,.35]:
            v=k.at(t);pred=48*np.cosh(4*t);actual=v['commutator_squared']/lam
            v.update(coefficient=actual,first_order=float(pred),relative_coefficient_error=float(actual/pred-1),
                subtraction_condition=(1+abs(v['otoc']))/(1-v['otoc']),coefficient_truncation_bound=2*cut['otoc_truncation_bound']/lam)
            points.append(v)
        refined=Kernel(lam,cut['L']+16).at(.35);diff=abs(refined['otoc']-points[-1]['otoc']);assert diff<1e-10
        rows.append(dict(lambda_value=lam,cutoff=cut,points=points,refined_cutoff=cut['L']+16,refinement_absolute_difference=diff,
            refinement_relative_difference=diff/abs(refined['otoc'])))
        print(lam,cut['L'],points[-1]['relative_coefficient_error'],round(time.monotonic()-start,1),flush=True)
    checks=[];mp.mp.dps=50
    for lam in ['.25','.1']:
        k=Kernel(float(lam),18)
        for t in ['.2','.7']:
            val,cond=nonsymmetric(mp.mpf(lam),mp.mpf(t),18);v=k.at(float(t));diff=abs(mp.mpf(v['otoc'])-val)
            assert diff<mp.mpf('1e-12') and abs(val.imag)<mp.mpf('1e-45')
            checks.append(dict(lambda_value=lam,time=t,L=18,dps=50,nonsymmetric_value=str(val),
                absolute_difference=str(diff),relative_difference=str(diff/abs(val)),sum_condition=str(cond)))
    extrap=[]
    for j,t in enumerate([0.,.1,.2,.35]):
        pred=rows[-1]['points'][j]['first_order'];value=2*rows[-1]['points'][j]['coefficient']-rows[-2]['points'][j]['coefficient']
        extrap.append(dict(time=t,value=value,first_order=pred,relative_error=value/pred-1))
    fig,ax=plt.subplots(figsize=(6.8,3.5));tt=np.linspace(0,.36,120)
    ax.plot(tt,48*np.cosh(4*tt),'k-',label='fixed-time leading coefficient')
    for r in rows:ax.plot([p['time'] for p in r['points']],[p['coefficient'] for p in r['points']],'o--',label=f"lambda = {r['lambda_value']}")
    ax.set(xlabel=r'$\mathcal{J}_s t$',ylabel=r'$D_{23}(t)/\lambda$');ax.legend(fontsize=8);fig.tight_layout();fig.savefig(OUT/'mixed-growth.png',dpi=180);plt.close(fig)
    inputs=['code/run_commutator_primary.py','code/mixed_cubic_otoc.py','code/singlet_otoc.py']
    data=dict(status='passed',scope='Numerical representation checks passed; finite-lambda errors measured, not asserted zero',
        conventions='Js=1, infinite temperature normalized trace, DS before fixed-time small lambda',records=rows,
        nonsymmetric_checks=checks,richardson=extrap,
        error_scope='Analytic infinite-cutoff Taylor bound; 50-digit original transfer checks at same finite L; no interval roundoff or uniform small-lambda remainder bound',
        environment=dict(python=platform.python_version(),numpy=np.__version__,mpmath=mp.__version__,matplotlib=matplotlib.__version__),
        input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs})
    target.write_text(json.dumps(data,indent=2)+'\n');print('saved',extrap,flush=True)
if __name__=='__main__':main()
