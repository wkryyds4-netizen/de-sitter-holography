"""Cross-representation diagnostics; preserve failed refinement and prototype bias."""
import hashlib
import json
import platform
import time
from pathlib import Path
import mpmath as mp
from ds3_racah_collision import physical_resolved,omega_collision,omega_collision_derivative,collision_data
from ds3_racah_merged import physical_merged,omega_merged
from global_racah import omega_quad
from run_ds3_collision import clean

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/stage4/ds3-collision'


def main():
    dest=BASE/'validation.json'
    if dest.exists():raise RuntimeError('Preserve prior run')
    data=dict(status='running',checks=[],heavy=[],versions=dict(python=platform.python_version(),mpmath=mp.__version__))
    def save():dest.write_text(json.dumps(clean(data),indent=2),encoding='utf-8')
    def check(name,value,reference,tol):
        error=abs(value-reference);rel=error/max(abs(reference),mp.mpf('1e-60'))
        data['checks'].append(dict(name=name,value=value,reference=reference,absolute_error=error,
                                   relative_error=rel,tolerance=mp.mpf(tol),passed=rel<mp.mpf(tol)))
        save();print(name,mp.nstr(rel,5),flush=True)
    mp.mp.dps=70
    # All three surfaces are imposed algebraically, not by decimal perturbation.
    for b0,d0,surface in [('.5',None,'a+c=b+d'),('.5','1.8','d=a+c+b'),('1.8','.5','b=a+c+d')]:
        a,c,p=map(mp.mpf,['.7','.6','.6']);b=mp.mpf(b0)
        d=a+c-b if d0 is None else (a+c+b if surface.startswith('d=') else mp.mpf(d0))
        if surface.startswith('b='):b=a+c+d
        v,detail=physical_merged(b,d,[a,p,c,p],128,28)
        ref,_=physical_resolved(b,d,[a,p,c,p],method='derivative')
        check(surface,v,ref,'1e-40')
        rev,_=physical_merged(d,b,[a,p,c,p],128,28)
        check(surface+'_Hermitian',v,mp.conj(rev),'1e-40')
        if surface=='a+c=b+d':
            hs=[(1+1j*x)/2 for x in [a,p,c,p]];H=(1-1j*b)/2;hp=(1+1j*d)/2
            old,oldmeta=omega_collision(H,hs,128,28);correct,newmeta=omega_merged(H,hp,hs,128,28)
            derivative=omega_collision_derivative(H,hs)
            data['prototype_bias']=dict(old=old,corrected=correct,reference=derivative,
                                       old_error=abs(old-derivative),corrected_error=abs(correct-derivative),
                                       cause='Python float division in digamma Bernoulli coefficient; corrected module uses mp arithmetic')
            # Independent direct Euler integral; endpoint numerical loss is recorded.
            with mp.workdps(45):euler=omega_quad(H,hp,hs)
            check('independent_Euler',euler,derivative,'1e-20')
            for off in ['1e-8','1e-25']:
                dd=d+mp.mpf(off);v,_=physical_merged(b,dd,[a,p,c,p],128,28)
                ref,_=physical_resolved(b,dd,[a,p,c,p])
                check('near_surface_'+off,v,ref,'1e-36')
    study=json.loads((BASE/'study.json').read_text())
    for mass in [8,16,32]:
        with mp.workdps(110 if mass==32 else 80):
            c=mp.mpf(mass);a=c+mp.mpf('.1');b=c+mp.mpf('.3');d=a+c-b;p=mp.mpf('.6')
            if mass<32:
                row=next(r for r in study['rows'] if r['M']==mass and r['offset']=='0')
                reference=mp.mpc(*row['runs'][-1]['value']);reference_source='study.json precision80; prototype bias below displayed scale here'
            else:
                reference,_=physical_resolved(b,d,[a,p,c,p],1024,36)
                reference_source='independent LPRS hyper for regular Omegas; exact-surface series for colliding Omega, dps110 N1024 K36'
            result=dict(M=mass,reference=reference,reference_source=reference_source,attempts=[])
            settings=[(max(128,16*mass),24),(max(256,32*mass),32)]
            if mass==32:settings += [(2048,40),(4096,48),(8192,56)]
            for N,K in settings:
                t=time.time();v,detail=physical_merged(b,d,[a,p,c,p],N,K)
                err=abs(v-reference)/abs(reference)
                result['attempts'].append(dict(N=N,K=K,value=v,relative_error=err,seconds=time.time()-t,
                                               shadow_condition=detail['shadow_condition'],
                                               omega_conditions=[r['summation_condition'] for r in detail['omega']]))
                print('heavy',mass,N,K,'relative',mp.nstr(err,6),flush=True)
            data['heavy'].append(result);save()
    data['status']='completed'
    paths=['code/validate_ds3_collision.py','code/ds3_racah_collision.py','code/ds3_racah_merged.py',
           'code/global_racah.py','code/run_ds3_collision.py','results/stage4/ds3-collision/study.json']
    data['hashes']={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}
    save()


if __name__=='__main__':main()
