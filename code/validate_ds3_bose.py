"""Independent finite Fock checks and continuous spherical derivative diagnostics."""
import hashlib
import itertools
import json
import platform
from pathlib import Path
import numpy as np
import mpmath as mp
from ds3_racah_fast import physical
from ds3_crossing_dictionary import physical_scalar

ROOT=Path(__file__).resolve().parents[1]


def main():
    out=ROOT/'results/stage4/ds3-bose/validation.json'
    if out.exists():raise RuntimeError('Preserve prior diagnostics')
    rng=np.random.default_rng(20261031);checks=[]
    def check(name,residual,tolerance):
        checks.append(dict(name=name,residual=float(residual),tolerance=float(tolerance),passed=bool(residual<=tolerance)))
    m=4;nobs=3
    Ms=rng.normal(size=(m,nobs,nobs))+1j*rng.normal(size=(m,nobs,nobs))
    f=rng.normal(size=nobs)+1j*rng.normal(size=nobs);f/=np.linalg.norm(f)
    h=np.array([.8,1.3,2.1]);basis=[v for v in itertools.product(range(3),repeat=m) if sum(v)<=2]
    ix={v:i for i,v in enumerate(basis)};creators=[]
    for j in range(m):
        A=np.zeros((len(basis),len(basis)))
        for k,v in enumerate(basis):
            w=list(v);w[j]+=1;w=tuple(w)
            if w in ix:A[ix[w],k]=np.sqrt(v[j]+1)
        creators.append(A)
    vac=np.zeros(len(basis));vac[ix[(0,)*m]]=1
    def mt(t):return Ms*np.exp(-1j*t*(h[:,None]-h[None,:]))[None,:,:]
    def labeled(t2,t1):return np.einsum('jab,kbc,c->ajk',mt(t2),mt(t1),f)
    def fock(t2,t1):
        A2=sum(np.kron(M,A) for M,A in zip(mt(t2),creators))
        A1=sum(np.kron(M,A) for M,A in zip(mt(t1),creators))
        return (A2@A1@np.kron(f,vac)).reshape(nobs,len(basis))
    def tensor_to_fock(B):
        result=np.zeros((nobs,len(basis)),complex)
        for j in range(m):
            for k in range(j,m):
                occupation=[0]*m;occupation[j]+=1;occupation[k]+=1
                result[:,ix[tuple(occupation)]]=B[:,j,k]*(1 if j==k else np.sqrt(2))
        return result
    for t2,t1 in [(0,0),(.7,0),(0,.7),(1.1,-.3)]:
        C=labeled(t2,t1);F=C.swapaxes(1,2);B=(C+F)/np.sqrt(2)
        actual=fock(t2,t1)
        check(f'fock_matrix_{t2}_{t1}',np.linalg.norm(actual-tensor_to_fock(B)),1e-12)
        check(f'direct_exchange_norm_{t2}_{t1}',abs(np.vdot(actual,actual)-np.vdot(C,C)-np.vdot(C,F)),1e-11)
    C=labeled(.7,0);D=labeled(0,.7);delta=C-D
    db=np.vdot(delta,delta)+np.vdot(delta,delta.swapaxes(1,2))
    check('order_difference_direct_fock',abs(db-np.linalg.norm(fock(.7,0)-fock(0,.7))**2),1e-11)
    check('equal_time_zero',np.linalg.norm(fock(0,0)-fock(0,0)),0)
    # Independent finite group average: permutation commutes with U_obs x U x U.
    charges=np.array([0,1,2,1]);qobs=np.array([0,1,2])
    total=qobs[:,None,None]+charges[None,:,None]+charges[None,None,:]
    eta=sum(np.exp(2j*np.pi*k*total/3) for k in range(3))/3
    eta=np.real_if_close(eta)
    etaB=eta*(C+C.swapaxes(1,2))/np.sqrt(2)
    rhs=np.vdot(C,eta*C)+np.vdot(C,eta*C.swapaxes(1,2))
    check('finite_group_descent',abs(np.vdot(etaB,etaB)-rhs),1e-11)
    check('finite_group_flip_commutation',np.max(abs(eta-eta.swapaxes(1,2))),1e-14)
    # Scalar observer (commuting coefficients) benchmark: permanent of Gram matrix.
    u=rng.normal(size=m)+1j*rng.normal(size=m);v=rng.normal(size=m)+1j*rng.normal(size=m)
    B=(np.outer(u,v)+np.outer(v,u))/np.sqrt(2)
    check('free_two_particle_permanent',abs(np.vdot(B,B)-(np.vdot(u,u)*np.vdot(v,v)+abs(np.vdot(u,v))**2)),1e-11)
    Bsame=np.sqrt(2)*np.outer(u,u)
    check('free_equal_vector_factor_two',abs(np.vdot(Bsame,Bsame)-2*np.vdot(u,u)**2),1e-11)
    # High-precision integral: positivity on a=c follows by pairing b=c+-x.
    mp.mp.dps=45;p=mp.mpf('.6');s=mp.mpf('.2')
    hh=lambda x:mp.sqrt(1+x*x)
    rr=lambda x:1/(mp.cosh(mp.pi*x)+mp.cosh(mp.pi*p))
    H=lambda a,b:rr(a-b)-rr(a+b)
    spherical=[]
    for c in map(mp.mpf,['.4','.6','.8']):
        integrand=lambda b:2*H(c,b)**2*(hh(b)-hh(c))
        direct=mp.quad(integrand,[s,c,2*c-s,2,5,mp.inf])
        paired=lambda x:2*(H(c,c+x)**2*(hh(c+x)-hh(c))-H(c,c-x)**2*(hh(c)-hh(c-x)))
        positive=mp.quad(paired,[0,c-s])+mp.quad(integrand,[2*c-s,2,5,mp.inf])
        check('diagonal_pairing_'+str(c),abs(direct-positive),mp.mpf('1e-40'))
        check('diagonal_strict_positive_'+str(c),0 if positive>0 else 1,0)
        spherical.append(dict(c=str(c),integral=str(direct),paired_integral=str(positive)))
    samples=[dict(a=.7,b=.5,c=.6,d=.8,kind='regular'),dict(a=.6,b=.7,c=.6,d=.7,kind='equal_labels')]
    with np.load(ROOT/'results/stage4/ds3-bose/pilot-n4-L8.npz') as z:
        for kind,index in [('max_condition',np.unravel_index(z['condition'].argmax(),z['condition'].shape))]:
            i,j=index;ia,ib,id=np.unravel_index(j,(len(z['x']),)*3)
            samples.append(dict(a=float(z['x'][ia]),b=float(z['x'][ib]),c=float(z['c'][i]),d=float(z['x'][id]),kind=kind))
    refs=[]
    for row in samples:
        a,b,c,d=[row[k] for k in ['a','b','c','d']]
        val,condition=physical(b,d,[a,.6,c,.6]);rev,_=physical(d,b,[a,.6,c,.6])
        ref=physical_scalar(mp.mpf(b),mp.mpf(d),list(map(mp.mpf,[a,.6,c,.6])))
        err=abs(mp.mpc(complex(val))-ref)
        refs.append(dict(**row,value=[float(val.real),float(val.imag)],reference=[str(ref.real),str(ref.imag)],
                         absolute_error=str(err),relative_error=str(err/abs(ref)),condition=float(condition),
                         hermiticity_residual=float(abs(val-np.conj(rev)))))
        print(row['kind'],'relative',str(err/abs(ref)),flush=True)
        if row['kind']!='max_condition':check('physical_flip_hermiticity_'+row['kind'],abs(val-np.conj(rev)),1e-10)
    paths=['code/validate_ds3_bose.py','code/ds3_racah_fast.py','code/ds3_crossing_dictionary.py',
           'code/global_racah.py','results/stage4/ds3-bose/pilot-n4-L8.npz']
    result=dict(status='passed' if all(c['passed'] for c in checks) else 'failed',checks=checks,
                finite_model='3 observer levels, 4 matter modes, 0/1/2 occupations; C3 group test is not dS group averaging',
                seed=20261031,versions=dict(python=platform.python_version(),numpy=np.__version__,mpmath=mp.__version__),
                dps=45,spherical_derivatives=spherical,kernel_references=refs,
                hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths})
    out.write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(result['status'],len(checks),flush=True)


if __name__=='__main__':main()
