"""Independent boundary witness and positive-tail checks, with actual raw values."""
import hashlib
import json
import platform
from pathlib import Path
import mpmath as mp
import numpy as np
from scipy.integrate import quad
from scipy.special import sici,roots_legendre
from ds3_group_pairing import emission000
from ds3_crossed_packet import positive_norms,probability
from ds3_packet_tail import bounds
from ds3_cutoff_witness import witness

ROOT=Path(__file__).resolve().parents[1]


def run():
    mp.mp.dps=45;checks=[];emissions=[];commutators=[];tails=[]
    def check(name,observed,expected,tol):
        observed,expected=mp.mpc(observed),mp.mpc(expected)
        err=abs(observed-expected)/max(abs(expected),mp.mpf('1e-30'))
        checks.append(dict(name=name,observed=str(observed),expected=str(expected),
            error=str(err),tolerance=str(tol),passed=bool(err<mp.mpf(tol))))
    def H(a,b,p):return 1/(mp.cosh(mp.pi*(a-b))+mp.cosh(mp.pi*p))-1/(mp.cosh(mp.pi*(a+b))+mp.cosh(mp.pi*p))
    def C(p):return mp.sqrt(mp.sinh(mp.pi*p)/p)/(4*mp.sqrt(mp.pi))
    for aa,bb,pp in [('.4','.65','.6'),('1.3','.7','.9'),('2.1','.25','1.2')]:
        a,b,p=map(mp.mpf,[aa,bb,pp])
        old=emission000(a,b,p);new=C(p)*mp.sqrt(mp.cosh(mp.pi*b)/mp.cosh(mp.pi*a))*H(a,b,p)
        emissions.append(dict(a=str(a),b=str(b),p=str(p),old=str(old),new=str(new)))
        check('emission '+aa,old,new,'1e-38')
    p,q=mp.mpf('.6'),mp.mpf('.9');s=mp.mpf('.2')
    for aa,cc in [('1.2','.6'),('1.1','.8'),('.5','.7')]:
        a,c=map(mp.mpf,[aa,cc])
        integrand=lambda b:H(a,b,q)*H(b,c,p)-H(a,b,p)*H(b,c,q)
        direct=mp.quad(integrand,[s,1,2,4,8,14,20])
        missing=-mp.quad(integrand,[0,s])
        tail=4/mp.pi*mp.exp(-mp.pi*(40-a-c))
        commutators.append(dict(a=str(a),c=str(c),s=str(s),direct=str(direct),missing=str(missing),tail=str(tail)))
        check('cutoff commutator '+aa,direct,missing,'1e-34')
    norms=[positive_norms(order=n,nc=16,L=20) for n in [20,32]]
    check('positive ns refinement',norms[0]['ns'],norms[1]['ns'],'1e-12')
    check('positive nt refinement',norms[0]['nt'],norms[1]['nt'],'1e-12')
    check('initial normalization',norms[1]['initial_norm'],1,'1e-13')
    for L in [8,10,12,16,20]:
        r=bounds(L);tails.append({k:str(v) for k,v in r.items()})
        # Direct elementary integration of the derived two upper envelopes.
        pi=mp.pi;gmin=r['gmin'];C0=r['C0'];CL=mp.exp(pi/(mp.sqrt(1+L*L)+L))
        ta=CL/(16*pi*pi*gmin)*mp.quad(lambda a:a*(a-mp.mpf('.4')+1/pi)*mp.exp(-pi*(a-mp.mpf('.8'))),[L,mp.inf])
        tb=1/(16*pi*pi*gmin)*mp.quad(lambda b:(2*(b+C0)/pi+2/pi**2)*mp.exp(-pi*(b-mp.mpf('.8'))),[L,mp.inf])
        check('tail a '+str(L),ta,r['tail_a'],'1e-30')
        check('tail b '+str(L),tb,r['tail_b'],'1e-30')
    witnesses=[witness(n) for n in [16,24,40]]
    for r in witnesses[:2]:
        check('witness overlap '+str(r['order']),r['overlap'],witnesses[-1]['overlap'],'1e-12')
        check('witness norm '+str(r['order']),r['test_norm_squared'],witnesses[-1]['test_norm_squared'],'1e-12')
    # Independent direct radial group average for the test vector norm.
    xs,ws=roots_legendre(48);cs=.4+.4*(xs+1)/2;cw=ws*.2
    f2=probability(cs)/(2*cs*cs/np.pi)
    def radial(r):
        if r==0:return 0.
        ia=(sici(1.4*r)[0]-sici(1.1*r)[0])/np.sinh(r)
        ic=np.sum(cw*f2*np.sin(cs*r)/(cs*np.sinh(r)))
        return np.sin(.6*r)*np.sin(.9*r)/(.6*.9)*ia*ic
    radial_value=sum(quad(radial,l,l+2,epsabs=1e-13,epsrel=1e-12)[0] for l in range(0,40,2))
    # |sin|<=1, sinh^-2 r <=4 e^-2r/(1-e^-80)^2 for r>=40.
    radial_tail=2*np.exp(-80)/(1-np.exp(-80))**2*np.log(1.4/1.1)*np.sum(cw*f2/cs)/(.6*.9)
    check('independent radial witness norm',radial_value,witnesses[-1]['test_norm_squared'],'1e-11')
    paths=['code/validate_ds3_packet_bounds.py','code/ds3_packet_tail.py',
           'code/ds3_cutoff_witness.py','code/ds3_crossed_packet.py','code/ds3_group_pairing.py']
    data=dict(status='passed' if all(c['passed'] for c in checks) else 'failed',
              checks=checks,emissions=emissions,commutators=commutators,positive_norms=norms,
              tails=tails,witnesses=witnesses,radial=dict(value=radial_value,tail_bound=float(radial_tail),L=40,order=48),
              versions=dict(python=platform.python_version(),mpmath=mp.__version__,numpy=np.__version__),
              hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths})
    dest=ROOT/'results/stage4/ds3-crossed-packet/bounds-validation.json'
    dest.write_text(json.dumps(data,indent=2),encoding='utf-8')
    print(data['status'],len(checks),[c for c in checks if not c['passed']],flush=True)


if __name__=='__main__':run()
