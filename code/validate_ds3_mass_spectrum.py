"""Independent spectral completeness, wavepacket and normalization diagnostics."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,platform
import mpmath as mp
import numpy as np
from ds3_mass_spectrum import *
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/stage4/ds3-mass-spectrum'

def encode(x):
    if isinstance(x,dict):return {k:encode(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [encode(v) for v in x]
    if isinstance(x,(mp.mpf,mp.mpc)):return mp.nstr(x,60)
    if isinstance(x,np.generic):return x.item()
    return x

def packet(T,nmass,ntime):
    # A smooth-enough compact mass packet, not a delta-normalization point sample.
    a,b=.6,1.8
    x,w=np.polynomial.legendre.leggauss(nmass)
    nu=(a+b)/2+(b-a)/2*x;w=w*(b-a)/2
    h=(nu-a)**4*(b-nu)**4
    h=h/np.sqrt(np.dot(w,h*h))
    q=np.exp(-np.pi*nu);den=np.sqrt(1+q*q)
    xm,wm=np.polynomial.legendre.leggauss(ntime)
    t=T*xm;wt=T*wm
    e=np.exp(-1j*np.outer(t,nu))
    B=(e@(w*h/den)-e.conj()@(w*h*q/den))/np.sqrt(2*np.pi)
    A=(e.conj()@(w*h/den)-e@(w*h*q/den))/np.sqrt(2*np.pi)
    O=(e@(w*h*q/den)+e.conj()@(w*h/den))/np.sqrt(2*np.pi)
    return dict(T=T,nmass=nmass,ntime=ntime,norm=float(np.dot(wt,abs(B)**2)),
        antipodal_overlap=complex(np.dot(wt,B.conj()*A)),
        orthogonal_overlap=complex(np.dot(wt,B.conj()*O)),
        expected_antipodal=float(-np.dot(w,h*h/np.cosh(np.pi*nu))))

def run():
    mp.mp.dps=50;checks=[];modes=[];bounds=[];grams=[];parseval=[];packets=[];legs=[]
    def check(name,value,expected,tol='1e-38',scale=None):
        value=mp.mpc(value);expected=mp.mpc(expected)
        scale=max(abs(expected),mp.mpf(1)) if scale is None else mp.mpf(scale)
        error=abs(value-expected)/scale
        checks.append(dict(name=name,observed=value,expected=expected,scale=scale,
            scaled_error=error,tolerance=tol,passed=bool(error<mp.mpf(tol))))
    swap=mp.matrix([[0,1],[1,0]])
    for ns in ['.03','.4','1.1','3.2']:
        nu=mp.mpf(ns);q=mp.exp(-mp.pi*nu);N=1/mp.sqrt(2*nu*(1-q*q))
        R=mp.sqrt(2*mp.pi)*N*mp.matrix([[1,-q],[-q,1]])
        G=gram(nu);P=projector(nu);v=bd_vector(nu);o=orth_vector(nu)
        check('Gram_'+ns,mp.norm(R.T*R-G),0)
        check('Gram_det_'+ns,mp.det(G),(mp.pi/nu)**2)
        check('Gram_inverse_'+ns,mp.norm(R*(G**-1)*R.T-mp.eye(2)),0)
        check('P_square_'+ns,mp.norm(P*P-P),0)
        check('P_orthogonal_'+ns,(v.T*o)[0],0)
        check('P_complete_'+ns,mp.norm(P+o*o.T-mp.eye(2)),0)
        check('antipodal_overlap_'+ns,(v.T*swap*v)[0],-mp.sech(mp.pi*nu))
        # Explicitly verify the antipodal line is NOT the orthogonal complement.
        defect=mp.norm(P+swap*P*swap-mp.eye(2))
        check('antipodal_defect_'+ns,defect,mp.sqrt(2)*mp.sech(mp.pi*nu))
        grams.append(dict(nu=nu,a=G[0,0],b=G[0,1],determinant=mp.det(G),
            eigenvalues=[G[0,0]+G[0,1],G[0,0]-G[0,1]],mu=mu(nu),
            overlap=(v.T*swap*v)[0],projector=[list(P[i,:]) for i in range(2)],defect=defect))
        for l in [0,1,3]:
            for ts in ['-.7','.2']:
                t=mp.mpf(ts);em=continuum(nu,l,t,-1);ep=continuum(nu,l,t,1)
                bd=mp.cosh(t)*mode(nu,l,t)
                check(f'BD_expand_{ns}_{l}_{ts}',mp.sqrt(mu(nu))*bd,v[0]*em+v[1]*ep)
                anti=(-1)**l*mp.cosh(t)*mode(nu,l,-t)
                check(f'antipode_expand_{ns}_{l}_{ts}',mp.sqrt(mu(nu))*anti,v[1]*em+v[0]*ep)
                modes.append(dict(nu=nu,l=l,t=t,bd=bd,antipode=anti,e_minus=em,e_plus=ep))
    print('Gram and channel tests complete',flush=True)
    for l in [1,2,3,4]:
        for k in range(1,l+1):
            f=lambda t:bound(l,k,t)
            norm=mp.quad(lambda t:abs(f(t))**2,[-mp.inf,-2,0,2,mp.inf])
            check(f'bound_norm_{l}_{k}',norm,1)
            t=mp.mpf('.37')
            residual=-mp.diff(f,t,2)-l*(l+1)*mp.sech(t)**2*f(t)+k*k*f(t)
            check(f'bound_ODE_{l}_{k}',residual,0)
            # Continuum-bound orthogonality checked by independent real-time quadrature.
            nu=mp.mpf('.8')
            overlap=mp.quad(lambda t:mp.conj(continuum(nu,l,t))*f(t),[-mp.inf,-4,0,4,mp.inf])
            check(f'bound_continuum_{l}_{k}',overlap,0,'1e-30')
            bounds.append(dict(l=l,k=k,norm=norm,ode_residual=residual,overlap=overlap))
        for k in range(1,l):
            overlap=mp.quad(lambda t:bound(l,k,t)*bound(l,l,t),[-mp.inf,-2,0,2,mp.inf])
            check(f'bound_cross_{l}_{k}',overlap,0)
    print('bound sector complete',flush=True)
    for ps,l in [('.5',0),('.5',1),('1.3',2),('2.2',3)]:
        p=mp.mpf(ps)
        overlaps=[sech_bound_overlap(p,l,k) for k in range(1,l+1)]
        bound_weight=mp.fsum(abs(a)**2 for a in overlaps)
        continuous=2*mp.quad(lambda nu:abs(sech_transform(p,l,nu))**2,[0,1,3,8,mp.inf])
        expected=mp.beta(p,mp.mpf('.5'))
        check(f'Parseval_{ps}_{l}',continuous+bound_weight,expected,'1e-32')
        # Independent transform at one nonzero frequency, not just a norm identity.
        nu=mp.mpf('.9')
        direct=mp.quad(lambda t:mp.conj(continuum(nu,l,t))*mp.sech(t)**p,
                       [-mp.inf,-6,0,6,mp.inf])
        analytic=sech_transform(p,l,nu)
        check(f'Fourier_direct_{ps}_{l}',direct,analytic,'1e-28')
        parseval.append(dict(p=p,l=l,continuum=continuous,bound_weight=bound_weight,
            bound_overlaps=overlaps,total=continuous+bound_weight,expected=expected,
            transform_direct=direct,transform_beta=analytic))
    print('independent Parseval checks complete',flush=True)
    for T in [20,40,80]:
        row=packet(T,140,360);ref=packet(T,220,540)
        for key in ['norm','antipodal_overlap','orthogonal_overlap']:
            check(f'packet_refine_{T}_{key}',row[key],ref[key],'2e-12')
        check(f'packet_norm_{T}',ref['norm'],1,'2e-7')
        check(f'packet_antipodal_{T}',ref['antipodal_overlap'],ref['expected_antipodal'],'2e-7')
        check(f'packet_orthogonal_{T}',ref['orthogonal_overlap'],0,'2e-7')
        for record in [row,ref]:
            for key in ['antipodal_overlap','orthogonal_overlap']:
                record[key]=mp.mpc(record[key])
        packets.append(dict(coarse=row,refined=ref))
    for ss in [['.7','1.1','.6'],['3.2','2.7','.5']]:
        nus=list(map(mp.mpf,ss))
        for cc in [(True,False,False),(True,False,True)]:
            c=observer_coefficient(nus,cc)
            old=reduced_coefficient(nus,cc)
            check(f'observer_leg_{ss}_{cc}',abs(c)**2/abs(old)**2,mu(nus[0])*mu(nus[1]))
            legs.append(dict(nus=nus,conjugates=cc,old_coefficient=old,coefficient=c,
                             norm_ratio=abs(c)**2/abs(old)**2,expected=mu(nus[0])*mu(nus[1])))
    paths=['code/ds3_mass_spectrum.py','code/validate_ds3_mass_spectrum.py','code/ds3_local_vertex.py',
        'sources/lagogiannis-maloney-wang-1106.2846v1.html','sources/marolf-morrison-0810.5163v1.html',
        'sources/kolchmeyer-liu-2411.08090v1.html','notes/stage4-ds3-local-vertex.md']
    data=dict(status='passed' if all(r['passed'] for r in checks) else 'failed',
        created_at_utc=datetime.now(timezone.utc).isoformat(),python=platform.python_version(),
        mpmath=mp.__version__,numpy=np.__version__,dps=50,checks=checks,
        gram=grams,modes=modes,bounds=bounds,parseval=parseval,packets=packets,legs=legs,
        scope='Spacetime L2 spectral/projector model. Finite numerical quadratures are empirical, not interval certificates. No physical group averaging or finite G.',
        hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths})
    OUT.mkdir(parents=True,exist_ok=True);dest=OUT/'diagnostics.json'
    assert not dest.exists(),'Do not overwrite evidence'
    dest.write_text(json.dumps(encode(data),indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=data['status'],checks=len(checks),failed=[r['name'] for r in checks if not r['passed']])))
if __name__=='__main__':run()
