"""Explicit candidate mass average, not an asserted finite-kappa dictionary.

Uses CLP Z_obs and the stipulated semiclassical mu= kappa(1-alpha/2pi)
at finite kappa to define a comparison model. G stays fixed across mu.
"""
import hashlib,json
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.special import ndtri
from scipy.stats import qmc
from principal_eikonal import kernel,spectral_moment
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/stage4/mass-filter'

def density(mu,kappa):
    mu=np.asarray(mu);alpha=2*np.pi*(1-mu/kappa)
    valid=(mu>0)&(mu<kappa);safe=np.where(valid,alpha,1.)
    return np.where(valid,(.5*np.sinc(safe/(2*np.pi))*(-np.expm1(-kappa*safe)))**2,0.)

def run(power,seed,normalizer,base):
    M=1000.;kappa=4000.;delta=1/np.sqrt(M);sd=1/(np.sqrt(8)*delta);t=np.log(M)
    sampler=qmc.Sobol(9,scramble=True,seed=seed);sums=np.zeros(5,complex)
    for _ in range(2**(power-15)):
        a=np.clip(sampler.random(2**15),np.finfo(float).eps,1-np.finfo(float).eps)
        z=ndtri(a[:,:4]);X,Y,Z=(delta*z[:,:3]/np.sqrt(2)).T;r=-np.log(a[:,4:8])
        A=np.exp(Y)*np.cosh(r[:,0])+np.exp(-Y)*np.cosh(r[:,1]);B=np.exp(X)*np.cosh(r[:,2])+np.exp(-X)*np.cosh(r[:,3])
        w=np.exp(r.sum(axis=1))/(A*A*B*B)
        w*=np.cos(.6*r[:,0])*np.cos(.6*r[:,1])*np.cos(.9*r[:,2])*np.cos(.9*r[:,3])
        phi=np.pi*a[:,8];T=t+Z
        def integrand(mu):
            h=-np.cos((1-mu/kappa)*(np.pi-phi))/(2*mu*np.sinc(mu/kappa))
            return kernel(h*np.exp(T)/(A*B))+kernel(h*np.exp(-T)/(A*B),-1)-2
        c=integrand(M);num=0j;massweight=0.
        for sign in [-1,1]:
            mu=M+sign*sd*z[:,3];p=density(mu,kappa);safe=np.clip(mu,1e-8,kappa-1e-8)
            num=num+p*integrand(safe)/2;massweight=massweight+p/2
        packet=w*num/normalizer
        conditional=w*c
        # This difference control has exactly zero expectation if c is mass-independent.
        difference=w*(num-massweight*c)/normalizer
        sums+=np.array([packet.sum(),conditional.sum(),difference.sum(),massweight.sum(),np.abs(w+packet).sum()])
    q=sums/2**power;pair=lambda v:[float(v.real),float(v.imag)]
    return dict(power=power,seed=seed,packet=pair(1+q[0]/base),conditional=pair(1+q[1]/base),
                difference=pair(q[2]/base),sampled_mass_normalization=float(q[3].real),absolute_mass_average=float(q[4].real/base))

def main():
    if (OUT/'packet-runs.jsonl').exists():raise RuntimeError('Do not overwrite')
    M=1000.;kappa=4000.;delta=1/np.sqrt(M);sd=1/(np.sqrt(8)*delta)
    norm,err=quad(lambda z:float(density(M+sd*z,kappa))*np.exp(-z*z/2)/np.sqrt(2*np.pi),-12,12,epsabs=1e-14,epsrel=1e-13)
    base=float(spectral_moment(0,.6,delta,50)*spectral_moment(0,.9,delta,50))
    rows=[];summaries=[]
    with (OUT/'packet-runs.jsonl').open('x') as file:
        for power in [19,21]:
            group=[]
            for seed in range(8):
                row=run(power,17100+seed,norm,base);rows.append(row);group.append(row);file.write(json.dumps(row)+'\n');file.flush()
            s=dict(power=power)
            for key in ['packet','conditional','difference']:
                x=np.array([r[key] for r in group]);s[key]=x.mean(axis=0).tolist();s[key+'_se']=(x.std(axis=0,ddof=1)/np.sqrt(8)).tolist()
            summaries.append(s);print(json.dumps(s),flush=True)
    inputs=['code/packet_mass_filter.py','code/principal_eikonal.py']
    out=dict(status='completed; explicitly stipulated candidate mass average',parameters=dict(M=M,kappa=kappa,G=1/(4*kappa),delta=delta,u=0),
             normalizer=norm,normalizer_quad_error_estimate=err,normalizer_tail_bound=float(.25*__import__('scipy').special.erfc(12/np.sqrt(2))),
             base=base,rows=rows,summaries=summaries,
             caveat='Finite-kappa mass-angle map and insertion completion are not established; this is a model comparison, not a physical finite-G trace.',
             hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs})
    (OUT/'packet.json').write_text(json.dumps(out,indent=2))
if __name__=='__main__':main()
