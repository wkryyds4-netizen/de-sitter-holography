"""Independent spherical, group-average, local-gluing and two-point checks."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,platform
import mpmath as mp
import numpy as np
from ds3_group_pairing import *
from recoil_spectral import joint_kernel,spectrum,correction
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/stage4/ds3-group-pairing'

def encode(x):
    if isinstance(x,dict):return {k:encode(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [encode(v) for v in x]
    if isinstance(x,(mp.mpf,mp.mpc)):return mp.nstr(x,60)
    if isinstance(x,np.generic):return x.item()
    if isinstance(x,complex):return mp.nstr(mp.mpc(x),60)
    return x

def time_square(nu,s,t,order,L):
    xs,ws=np.polynomial.legendre.leggauss(order)
    a=s+(L-s)*(xs+1)/2;w=ws*(L-s)/2
    A=a[:,None];B=a[None,:]
    D=(np.cosh(np.pi*(A+B))+np.cosh(np.pi*nu))*(np.cosh(np.pi*(A-B))+np.cosh(np.pi*nu))/4
    Q=A*B*np.tanh(np.pi*A)*np.tanh(np.pi*B)/(8*np.pi**2*D)
    E=np.sqrt(1+A*A)-np.sqrt(1+B*B)
    return complex(np.sum(w[:,None]*w[None,:]*Q*np.exp(-1j*t*E)))/float(state_norm(mp.mpf(str(s))))

def run():
    mp.mp.dps=50;checks=[];spherical_rows=[];pair_rows=[];triple_rows=[];gluing=[];joint_rows=[];states=[];rates=[];times=[];means=[]
    def check(name,observed,expected,tol='1e-38',scale=None):
        observed=mp.mpc(observed);expected=mp.mpc(expected)
        scale=max(abs(expected),mp.mpf('1e-60')) if scale is None else mp.mpf(scale)
        error=abs(observed-expected)/scale
        checks.append(dict(name=name,observed=observed,expected=expected,scale=scale,scaled_error=error,tolerance=tol,passed=bool(error<mp.mpf(tol))))
    for ns in ['.4','1.1','2.7']:
        n=mp.mpf(ns)
        for rs in ['.2','1','3']:
            r=mp.mpf(rs)
            angular=mp.quad(lambda x:(mp.cosh(r)-mp.sinh(r)*x)**(-1-1j*n),[-1,0,1])/2
            exact=spherical(n,r)
            check(f'spherical_{ns}_{rs}',angular,exact)
            spherical_rows.append(dict(nu=n,r=r,angular=angular,exact=exact))
    a=mp.mpf('.4');b=mp.mpf('1.6');target=compact_pair_norm(a,b)
    for R in [40,80]:
        radial=mp.quad(lambda r:compact_sine(r,a,b)**2,list(range(0,R+1,2)))
        tail=compact_pair_tail(R,a,b)
        check(f'pair_packet_tail_{R}',radial-target,0,'1',scale=tail)
        pair_rows.append(dict(a=a,b=b,R=R,radial=radial,schur=target,tail_bound=tail,error=abs(radial-target)))
    print('spherical functions and two-observer packet complete',flush=True)
    samples=[['.7','1.1','.6'],['1.4','.9','1.2'],['3.2','.4','.5']]
    for ss in samples:
        nus=list(map(mp.mpf,ss));kg=triple_group(nus);kt=triple_tanh(nus)
        direct,tail=triple_direct(nus)
        check('triple_tanh_'+str(ss),kt,kg)
        check('triple_integral_'+str(ss),direct,kg,'1e-23')
        check('triple_tail_'+str(ss),direct-kg,0,'1',scale=tail+mp.mpf('1e-45'))
        triple_rows.append(dict(nus=nus,closed=kg,tanh=kt,direct=direct,tail_bound=tail))
        a,b,nu=nus
        for ls in [(0,0,0),(1,1,0),(2,1,1),(2,2,2),(3,2,1)]:
            vl,_=vertex(nus,ls,(True,False,True))
            vr,_=vertex([b,a,nu],[ls[1],ls[0],ls[2]],(True,False,True))
            left=(-1)**ls[1]*f_linear(b)*vl
            right=(-1)**ls[0]*f_linear(a)*vr
            check('gluing_'+str(ss)+'_'+str(ls),left,right,'1e-36')
            gluing.append(dict(nus=nus,ls=ls,left=left,right=right))
        E=emission000(a,b,nu)
        raw=abs(f_linear(b)*E)**2/kg;closed=joint(a,b,nu)
        check('joint_raw_'+str(ss),raw,closed)
        check('joint_symmetry_'+str(ss),closed,joint(b,a,nu))
        ratio=(mp.pi/8)*closed/joint_kernel(a,b,nu)
        check('recoil_n_'+str(ss),ratio,local_n(a,b))
        joint_rows.append(dict(nus=nus,emission=E,state_amplitude=f_linear(b),triple_norm=kg,
            raw_joint=raw,closed=closed,old_joint=joint_kernel(a,b,nu),ratio=ratio,n=local_n(a,b)))
    print('triple averaging, gluing and physical kernel complete',flush=True)
    for ss in ['.2','.8','1.5']:
        s=mp.mpf(ss);z=state_norm(s);direct=mp.quad(state_density,[s,s+1,s+4,mp.inf])
        check('state_norm_'+ss,direct,z)
        states.append(dict(cutoff=s,analytic=z,direct=direct))
    for bs in ['.4','1.2','4']:
        b=mp.mpf(bs);s=mp.mpf('.2');nu=mp.mpf('.6')
        value=rate(b,nu,s);bound=rate_bound(b,s)
        check('creation_rate_bound_'+bs,value,0,'1',scale=bound)
        rates.append(dict(b=b,nu=nu,cutoff=s,rate=value,bound=bound))
    for Ms in ['3','6','12']:
        M=mp.mpf(Ms);nu=mp.mpf('.6')
        for ws in ['0','.3','1.2']:
            w=mp.mpf(ws);a=M+w/2;b=M-w/2
            reduced=joint(a,b,nu)/state_density(M)
            expected=spectrum(w,nu)*mean_correction(M,w,nu)
            check(f'mean_spectrum_{Ms}_{ws}',reduced,expected)
            delta=mp.sqrt(1+a*a)-mp.sqrt(1+b*b)
            phase_bound=abs(w)/(2*min(a,b)**2)
            if w:check(f'proper_energy_bound_{Ms}_{ws}',delta-w,0,'1',scale=phase_bound)
            means.append(dict(M=M,w=w,nu=nu,reduced=reduced,expected=expected,
                correction=mean_correction(M,w,nu),proper_difference=delta,label_difference=w,energy_error_bound=phase_bound))
    for nu,s in [(.6,.2),(1.2,.8)]:
        L=s+8
        for t in [0,.5,2,5]:
            coarse=time_square(nu,s,t,80,L);refined=time_square(nu,s,t,120,L)
            check(f'time_refine_{nu}_{s}_{t}',coarse,refined,'2e-10')
            check(f'time_imaginary_{nu}_{s}_{t}',refined.imag,0,'1e-14',scale=1)
            times.append(dict(nu=nu,cutoff=s,L=L,t=t,orders=[80,120],coarse=coarse,refined=refined,
                analytic_outer_tail=square_tail(mp.mpf(str(L)),mp.mpf(str(s)))/state_norm(mp.mpf(str(s)))))
    paths=['code/ds3_group_pairing.py','code/validate_ds3_group_pairing.py','code/ds3_local_vertex.py',
        'code/ds3_mass_spectrum.py','code/recoil_spectral.py','code/recoil_poles.py',
        'sources/marolf-yu-2610.24073v1.html','sources/marolf-morrison-0810.5163v1.html',
        'sources/kolchmeyer-liu-2411.08090v1.html','sources/derkachov-spiridonov-1711.07073v2.html']
    data=dict(status='passed' if all(r['passed'] for r in checks) else 'failed',
        created_at_utc=datetime.now(timezone.utc).isoformat(),python=platform.python_version(),mpmath=mp.__version__,numpy=np.__version__,dps=50,
        checks=checks,spherical=spherical_rows,pair_packets=pair_rows,triples=triple_rows,gluing=gluing,joint=joint_rows,
        state_norms=states,rates=rates,means=means,time_integrals=times,
        scope='Declared G=0 local observer model, group-averaged vacuum and one-matter sectors; linear-gluing state. Not finite kappa, full Fock algebra, interval quadrature or proof of a complete holographic dual.',
        hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths})
    OUT.mkdir(parents=True,exist_ok=True);dest=OUT/'diagnostics.json';assert not dest.exists()
    dest.write_text(json.dumps(encode(data),indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=data['status'],checks=len(checks),failed=[r['name'] for r in checks if not r['passed']])))
if __name__=='__main__':run()
