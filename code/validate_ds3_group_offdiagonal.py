"""Independent non-spherical group averages, testing the rank-one normalization."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,platform
import mpmath as mp
from ds3_local_vertex import future_phase,vertex
from ds3_group_pairing import triple_group,spherical
ROOT=Path(__file__).resolve().parents[1]

def boost_l0(nu,l,r):
    if r==0:return mp.mpf(l==0)
    factor=mp.conj(future_phase(nu,l))*mp.sqrt(2*l+1)/2
    if r<mp.mpf('.3'):
        return factor*mp.quad(lambda x:mp.legendre(l,x)*(mp.cosh(r)-mp.sinh(r)*x)**(-1-1j*nu),[-1,0,1])
    coeff={0:[1],1:[0,1],2:[-mp.mpf('.5'),0,mp.mpf('1.5')]}[l]
    total=mp.mpc(0)
    for n,c in enumerate(coeff):
        total+=c/mp.sinh(r)**(n+1)*mp.fsum(mp.binomial(n,k)*mp.cosh(r)**(n-k)*(-1)**k*
            2*mp.sinh((k-1j*nu)*r)/(k-1j*nu) for k in range(n+1))
    return factor*total

def run():
    mp.mp.dps=45;rows=[];checks=[]
    def check(name,a,b,tol='1e-32',scale=None):
        scale=max(abs(b),mp.mpf('1e-60')) if scale is None else scale
        err=abs(a-b)/scale
        checks.append(dict(name=name,observed=mp.nstr(a,55),expected=mp.nstr(b,55),scale=mp.nstr(scale,55),
            scaled_error=mp.nstr(err,55),tolerance=tol,passed=bool(err<mp.mpf(tol))))
    for ss in [['.7','1.1','.6'],['1.4','.9','1.2']]:
        a,b,p=map(mp.mpf,ss)
        T=70
        v0,_=vertex([a,b,p],(0,0,0),(True,False,True))
        for l in [1,2]:
            v,_=vertex([a,b,p],(l,l,0),(True,False,True))
            invariant=mp.sqrt(2*l+1)*v/v0
            expected=triple_group([a,b,p])*invariant
            direct=mp.quad(lambda r:(-1)**l/mp.sqrt(2*l+1)*mp.sinh(r)**2*
                boost_l0(a,l,r)*boost_l0(b,l,r)*spherical(p,r),list(range(0,T+1,2)))
            tail=2*mp.sqrt(2*l+1)/p/(1-mp.exp(-2*T))*mp.exp(-T)*(T*T+2*T+2)
            check(f'group_offdiag_{ss}_{l}',direct,expected,'1e-22')
            check(f'group_offdiag_tail_{ss}_{l}',direct-expected,0,'1',scale=tail)
            for r in [mp.mpf('.8'),mp.mpf('2.3')]:
                angular=mp.conj(future_phase(a,l))*mp.sqrt(2*l+1)/2*mp.quad(
                    lambda x:mp.legendre(l,x)*(mp.cosh(r)-mp.sinh(r)*x)**(-1-1j*a),[-1,0,1])
                check(f'boost_angular_{ss}_{l}_{r}',boost_l0(a,l,r),angular)
            rows.append(dict(nus=ss,l=l,T=T,direct=mp.nstr(direct,55),expected=mp.nstr(expected,55),
                invariant=mp.nstr(invariant,55),spherical_group=mp.nstr(triple_group([a,b,p]),55),tail_bound=mp.nstr(tail,55)))
            print('non-spherical group average',ss,l,flush=True)
    paths=['code/validate_ds3_group_offdiagonal.py','code/ds3_group_pairing.py','code/ds3_local_vertex.py']
    data=dict(status='passed' if all(c['passed'] for c in checks) else 'failed',checks=checks,rows=rows,dps=45,
        created_at_utc=datetime.now(timezone.utc).isoformat(),python=platform.python_version(),mpmath=mp.__version__,
        scope='Four independent non-spherical K-singlet to spherical group integrals. Analytic outer tails; no interior interval certificate.',
        hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths})
    dest=ROOT/'results/stage4/ds3-group-pairing/offdiagonal.json';assert not dest.exists()
    dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=data['status'],checks=len(checks),failed=[c['name'] for c in checks if not c['passed']])))
if __name__=='__main__':run()
