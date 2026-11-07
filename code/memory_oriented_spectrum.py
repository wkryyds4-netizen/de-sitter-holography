"""Independent low moments and two representations of the mixed memory kernel."""
from pathlib import Path
from math import exp,sqrt,pi
import json,hashlib,platform
import numpy as np
import sympy as s
from singlet_coefficient import pairings
from memory_oriented_marks import graph_oriented
from quartic_orbits import cross,q
from singlet_otoc import taylor_tail

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-oriented-kernel'


def moments():
    def hm(n):
        return sum(q**sum(cross(a,b,c,d) for i,(a,b) in enumerate(pa) for c,d in pa[i+1:])
                   for pa in pairings(tuple(range(n))))
    rows=[]
    for n in [0,2,4]:
        Sn=hm(n+4)-2*hm(n+2)+hm(n)
        for a in range(n+1):
            b=n-a;P=0
            for pa in pairings(tuple(range(n))):
                cr=sum(cross(i,j,k,l) for v,(i,j) in enumerate(pa) for k,l in pa[v+1:])
                P+=q**(cr+2*sum(i<a<=j for i,j in pa))
            marks=((0,1),(a+2,a+3));sums=[0,0,0];count=0
            for pa in pairings(tuple(range(n+4))):
                if any(mark in pa for mark in marks):continue
                g=graph_oriented(n+4,pa,marks,(False,True));count+=1
                sums[0]+=q**g['crossings']*g['delta']
                sums[1]+=q**g['crossings']*g['alpha']*g['beta']
                sums[2]+=q**g['crossings']*g['alpha']
            u=(1-q)/(1+q)
            assert s.factor(sums[0]+4*q*P)==0
            assert s.factor(sums[1]-4*u*u*Sn-16*q*P/(1+q))==0
            assert s.factor(sums[2]-2*u*Sn)==0
            rows.append(dict(a=a,b=b,pairings=count,probe_polynomial=str(s.expand(P))))
    return rows


def curves():
    source=json.loads((OUT/'reverse-cross.json').read_text())
    ll,C=s.symbols('lam C',positive=True)
    v4=s.sympify(source['V4'],locals={'q':q,'lam':ll,'C':C})
    v2=s.sympify(source['V2'],locals={'q':q,'lam':ll,'C':C})
    times=np.array([0,.25,.5,1]);rows=[]
    for lam in [.25,.1,.025]:
        qq=exp(-4*lam);bound=2/sqrt(lam*(1-qq));n=max(8,int(bound)+2)
        while 4*taylor_tail(bound,n)>1e-13:n+=1
        weights=np.array([float(v.subs({ll:lam,q:qq,C:1/lam})) for v in [v4,v2]])
        def jac(dim):
            off=np.sqrt(-np.expm1(-4*lam*np.arange(1,dim))/(-np.expm1(-4*lam)*lam))
            en,U=np.linalg.eigh(np.diag(off,1)+np.diag(off,-1))
            vec=U@(U[0,:,None]*np.exp(-1j*en[:,None]*times[None,:]))
            return np.array([np.sum(qq**(k*np.arange(dim)[:,None])*abs(vec)**2,axis=0) for k in [4,2]])
        aa=jac(n);bb=jac(n+12);err=float(np.max(abs(aa-bb)));assert err<2e-11
        data=dict(lam=lam,Js=1,weights_V4_V2=weights.tolist(),V2_fraction=float(weights[1]/weights.sum()),
                  cutoff=n,refined_cutoff=n+12,analytic_normalized_truncation_bound=4*taylor_tail(bound,n),
                  refinement_max_absolute=err,kernels_q4_q2=bb.tolist(),
                  normalized_mixture=(weights@bb/weights.sum()).tolist())
        if lam in [.25,.1]:
            checks=[]
            for nodes in [128,192]:
                xx,wt=np.polynomial.legendre.leggauss(nodes);th=(xx+1)*pi/2;wt=wt*pi/2
                rho=np.ones(nodes)/(2*pi);k=0
                while qq**k/(1-qq)>1e-18:
                    rho*=abs(1-qq**k*np.exp(2j*th))**2*(1-qq**(k+1));k+=1
                energies=2*np.cos(th)/sqrt(lam*(1-qq));de=energies[:,None]-energies[None,:]
                values=[];masses=[]
                for degree in [4,2]:
                    r=qq**degree;num=1.;den=np.ones((nodes,nodes))
                    for j in range(k):
                        v=r*qq**j;num*=1-r*r*qq**j
                        den*=(1-2*v*np.cos(th[:,None]-th[None,:])+v*v)*(1-2*v*np.cos(th[:,None]+th[None,:])+v*v)
                    measure=(rho*wt)[:,None]*(rho*wt)[None,:]*num/den
                    values.append([float(np.sum(measure*np.cos(t*de))) for t in times]);masses.append(float(measure.sum()))
                error=float(np.max(abs(np.array(values)-bb)));assert error<2e-10
                checks.append(dict(nodes=nodes,kernels=values,masses=masses,max_absolute_error=error))
            data['independent_mehler']=checks
        rows.append(data)
    return dict(times=times.tolist(),rows=rows,small_lambda_limit=(1/np.cosh(2*times)**8).tolist(),
                error_scope='Analytic Jacobi truncation; floating arithmetic and Mehler quadrature checked by refinement, not interval certified.')


def main():
    target=OUT/'spectral-checks.json'
    if target.exists():raise FileExistsError(target)
    mm=moments();print('opposite orientation moments',len(mm),flush=True)
    cc=curves()
    names=['code/memory_oriented_spectrum.py','code/memory_oriented_marks.py','code/singlet_coefficient.py',
           'code/quartic_orbits.py','code/singlet_otoc.py',
           'results/susskind-memory-oriented-kernel/reverse-cross.json']
    target.write_text(json.dumps(dict(status='passed',mixed_moments=mm,kernel=cc,
        environment=dict(python=platform.python_version(),numpy=np.__version__,sympy=s.__version__),
        input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names}),indent=2)+'\n')
    print('passed; mixed kernels saved',flush=True)


if __name__=='__main__':main()
