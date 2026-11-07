"""Unequal-size endpoint cancellation and the existing r=q^4 probe kernel."""
from pathlib import Path
from itertools import combinations,permutations
from math import factorial,exp,sqrt,pi
import json,hashlib,platform
import numpy as np
import sympy as s
from quartic_orbits import ww,constant,den,cross,polynomial_dict,accumulate
from singlet_otoc import taylor_tail

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-quartic'
q=s.Symbol('q');lam=s.Symbol('lambda',positive=True)


def cross_certificate(size):
    ma=list(combinations(range(4),2));mb=list(combinations(range(size),2))
    if size==2:rw=[1+q];rc=-2*(1-q)
    else:
        d=(1+q*q)*(1+q+q*q)
        rw=[-q*(1+q+q*q),d,-q*(1+q+q*q)];rc=-2*(1-q)*(1+q**3)
    prod=[[polynomial_dict(a*b) for b in rw+[rc]] for a in ww+[constant]]
    records=[];fails=[];total=0
    for k in range(size+1):
        nex=4+size-2*k
        for aset in combinations(range(nex),4-k):
            colors=''.join('A' if j in aset else 'B' for j in range(nex))
            for split in range(nex+1):
                A=list(range(4));B=list(range(4+split,4+split+size))
                ext=list(range(4,4+split))+list(range(4+split+size,4+size+nex))
                ae=[e for e,c in zip(ext,colors) if c=='A'];be=[e for e,c in zip(ext,colors) if c=='B']
                seed=list(zip(A[:k],B[:k]))+list(zip(A[k:],ae))+list(zip(B[k:],be));words=set()
                for ap in permutations(A):
                    for bp in permutations(B):
                        pm=dict(zip(A,ap))|dict(zip(B,bp))
                        words.add(tuple(sorted(tuple(sorted((pm.get(a,a),pm.get(b,b)))) for a,b in seed)))
                assert len(words)==24*factorial(size)//factorial(k)
                total+=len(words);tt={};dd={}
                for pa in words:
                    pi={a:b for a,b in pa}|{b:a for a,b in pa}
                    cr=sum(cross(a,b,c,d) for j,(a,b) in enumerate(pa) for c,d in pa[j+1:])
                    aa=[2-4*cross(A[a],pi[A[a]],A[b],pi[A[b]]) for a,b in ma]+[1]
                    bb=[2-4*cross(B[a],pi[B[a]],B[b],pi[B[b]]) for a,b in mb]+[1]
                    t=tt.setdefault(cr,[[0]*len(bb) for _ in aa]);d=dd.setdefault(cr,[[0]*len(mb) for _ in ma])
                    for i in range(len(aa)):
                        for j in range(len(bb)):t[i][j]+=aa[i]*bb[j]
                    for i,(a,b) in enumerate(ma):
                        for j,(c,e) in enumerate(mb):d[i][j]+=4*int(pi[A[a]]==B[e] and pi[A[b]]==B[c])
                linear={};quadratic={}
                for cr,t in tt.items():
                    for i in range(len(t)):
                        for j in range(len(t[0])):accumulate(quadratic,prod[i][j],t[i][j],cr)
                for cr,t in dd.items():
                    for i in range(len(t)):
                        for j in range(len(t[0])):accumulate(linear,prod[i][j],t[i][j],cr)
                linear={j:v for j,v in linear.items() if v};quadratic={j:v for j,v in quadratic.items() if v}
                if linear or quadratic:fails.append(dict(k=k,colors=colors,split=split,linear=linear,quadratic=quadratic))
                records.append(dict(bridges=k,colors=colors,split=split,words=len(words),vanishes=not linear and not quadratic))
    print('cross size',size,'templates',len(records),'counterexamples',len(fails),flush=True)
    return dict(size=size,templates=len(records),words=total,records=records,counterexamples=fails)


def curves():
    records=[];times=[0,.25,.5,1]
    for ll in [.25,.1,.025,.01]:
        qq=exp(-4*ll);B=2/sqrt(ll*(1-qq));n=max(8,int(B)+2)
        while 4*taylor_tail(B,n)>1e-13:n+=1
        def jac(n):
            off=np.sqrt(-np.expm1(-4*ll*np.arange(1,n))/(-np.expm1(-4*ll)*ll))
            en,U=np.linalg.eigh(np.diag(off,1)+np.diag(off,-1))
            v=U@(U[0,:,None]*np.exp(-1j*en[:,None]*np.array(times)[None,:]))
            return np.sum(qq**(4*np.arange(n)[:,None])*abs(v)**2,axis=0)
        a=jac(n);b=jac(n+12);err=float(np.max(abs(a-b)));assert err<2e-11
        row=dict(lam=ll,cutoff=n,refined_cutoff=n+12,analytic_truncation_bound=4*taylor_tail(B,n),
                 refinement=err,values=b.tolist(),small_lambda=(1/np.cosh(2*np.array(times))**8).tolist())
        if ll in [.25,.1]:
            vals=[]
            for nodes in [128,192]:
                xx,wt=np.polynomial.legendre.leggauss(nodes);th=(xx+1)*pi/2;wt=wt*pi/2
                rho=np.ones(nodes)/(2*pi);num=1.;d=np.ones((nodes,nodes));k=0;r=qq**4
                while qq**k/(1-qq)>1e-18:
                    rho*=abs(1-qq**k*np.exp(2j*th))**2*(1-qq**(k+1))
                    v=r*qq**k;num*=1-r*r*qq**k
                    d*=(1-2*v*np.cos(th[:,None]-th[None,:])+v*v)*(1-2*v*np.cos(th[:,None]+th[None,:])+v*v);k+=1
                measure=(rho*wt)[:,None]*(rho*wt)[None,:]*num/d
                en=2*np.cos(th)/sqrt(ll*(1-qq));de=en[:,None]-en[None,:]
                val=[float(np.sum(measure*np.cos(t*de))) for t in times]
                assert max(abs(np.array(val)-b))<2e-10
                vals.append(dict(nodes=nodes,values=val,mass=float(measure.sum())))
            row['independent_mehler']=vals
        records.append(row)
    return dict(times=times,records=records)


def main():
    target=OUT/'cross-kernel.json'
    if target.exists():raise FileExistsError(target)
    crosses=[cross_certificate(j) for j in [2,3]]
    prior=json.loads((OUT/'quartic.json').read_text());vv=s.sympify(prior['residual_variance'],locals={'q':q,'lambda':lam})
    certificate=json.loads((OUT/'orbit-certificate.json').read_text());survivor=certificate['records'][-1]
    poly=lambda obj:sum(int(v)*q**int(k) for k,v in obj.items())
    surviving=(lam*poly(survivor['linear'])+lam**2*poly(survivor['quadratic']))/den**2
    assert s.factor(surviving-vv)==0
    norm_lead=s.limit(vv.subs(q,s.exp(-4*lam))/lam**3,lam,0,dir='+')
    assert norm_lead==s.Rational(768,5)
    data=curves();names=['code/quartic_cross_kernel.py','code/quartic_orbits.py','code/singlet_otoc.py',
                       'results/susskind-memory-quartic/quartic.json','results/susskind-memory-quartic/orbit-certificate.json']
    result=dict(status='passed' if not any(x['counterexamples'] for x in crosses) else 'counterexamples',
                cross_certificates=crosses,small_lambda_variance_coefficient=str(norm_lead),
                surviving_orbit_matches_gram=True,kernel=data,
                environment=dict(python=platform.python_version(),numpy=np.__version__,sympy=s.__version__),
                input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names})
    target.write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
