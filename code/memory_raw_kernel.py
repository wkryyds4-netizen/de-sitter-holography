"""Positive raw three-field matrix and weighted Jacobi/Mehler verification."""
from pathlib import Path
from math import exp,sqrt,pi
import json,hashlib,platform
import numpy as np
import sympy as s
from memory_raw_spectrum import coefficients,q,lam,C,x,y
from singlet_otoc import taylor_tail

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-raw-spectrum'


def assembly():
    w,f,r,z,e,a,h,bf,br,ss,uu=s.symbols('w f r z e a h bf br s u',real=True)
    M=s.Matrix([[1,0,0,0],[a*e,1,1,0],[h,bf*e,br*e,1]])
    mat=M*s.diag(w,f,r,z)*M.T
    det=w*((f+r)*z+e*e*f*r*(bf-br)**2)
    assert s.factor(mat.det()-det)==0
    density=(s.Matrix([[1,-s.I*ss,-ss**2/2]])*mat*s.Matrix([1,s.I*uu,-uu**2/2]))[0]
    even=s.expand((density+density.subs(e,-e))/2)
    expected=w+ss*uu*mat[1,1]-(ss**2+uu**2)*w*h/2+ss**2*uu**2*mat[2,2]/4
    assert s.expand(even-expected)==0
    return dict(matrix=[[str(v) for v in row] for row in mat.tolist()],determinant=str(det),
                full_polynomial_density=str(s.expand(density)),ordinary_polynomial_density=str(even))


def curves():
    saved=json.loads((OUT/'symbolic-spectrum.json').read_text());co=coefficients()
    parse=lambda z:s.sympify(z,locals={'q':q,'lam':lam,'C':C})
    norms=list(map(parse,saved['norm_parts']));curvs=list(map(parse,saved['curvature_parts']))
    prior=json.loads((ROOT/'results/susskind-memory-oriented-kernel/reverse-cross.json').read_text())
    times=np.array([0,.25,.5,.75,1]);records=[]
    for ll in [.25,.1,.025,.01]:
        llrat=s.Rational(str(ll));subs={lam:llrat,q:s.exp(-4*llrat),C:1/llrat}
        ev=lambda z:float(s.N(z.subs(subs),45))
        cc=1/ll;qq=exp(-4*ll);rad=2/sqrt(ll*(1-qq))
        bj,bd,bf,br,kap,aa=[ev(co[k]) for k in ['bJ','bD','bF','bR','kappa','a']]
        W=4/ll+16*qq/(1+qq);VF=64*qq*qq
        D=(1+qq*qq)*(1+qq+qq*qq)
        VR=cc**3*(ll*4*(1-qq*qq)*(1+qq*qq)*(1+qq+qq*qq)**2+
             ll**2*16*qq*qq*(1+qq)*(1+qq*qq)*(1+qq+qq*qq)*(2*qq*qq-qq+2))/D**2
        V4,V2=[ev(parse(prior[k])) for k in ['V4','V2']]
        norm=sum(ev(z) for z in norms);curvature=sum(ev(z) for z in curvs)
        K=bj*(x+y)**2/4+bd*(x-y)**2+cc*kap
        e2=(x+y)**2/4
        def pd(poly):return [(i,j,float(v)) for (i,j),v in s.Poly(s.expand(poly),x,y).terms()]
        pK,pK2,pe2,pomega2=pd(K),pd(K*K),pd(e2),pd((x-y)**2)
        kbound=(abs(bj)+4*abs(bd))*rad**2+abs(cc*kap)
        gain=V4+V2+(bf*bf*VF+br*br*VR)*rad**2+W*kbound*kbound
        cutoff=max(12,int(rad)+8)
        while 4*gain/V4*taylor_tail(rad,cutoff-4)>1e-12:cutoff+=1
        def jac(n):
            off=np.sqrt(-np.expm1(-4*ll*np.arange(1,n))/(-np.expm1(-4*ll)*ll))
            en,U=np.linalg.eigh(np.diag(off,1)+np.diag(off,-1))
            phase=np.exp(-1j*en[:,None]*times[None,:])
            vectors=[U@(U[0,:,None]*en[:,None]**j*phase) for j in range(5)]
            def moment(power,poly):
                rn=qq**(power*np.arange(n))
                terms=np.array([c*np.sum(rn[:,None]*vectors[i].conj()*vectors[j],axis=0) for i,j,c in poly])
                value=terms.sum(axis=0);assert max(abs(value.imag))<1e-8
                return value.real
            raw=V4*moment(4,[(0,0,1)])+V2*moment(2,[(0,0,1)])
            raw+=bf*bf*VF*moment(1,pe2)+br*br*VR*moment(3,pe2)+W*moment(2,pK2)
            yu=W*moment(2,pK)
            return np.array([raw/norm,yu/sqrt(W*norm)])
        base=jac(cutoff);fine=jac(cutoff+12);err=float(np.max(abs(base-fine)));assert err<2e-10
        assert abs(fine[0,0]-1)<2e-10
        row=dict(lam=ll,Js=1,cutoff=cutoff,refined_cutoff=cutoff+12,
                 normalized_raw_truncation_bound=4*gain/V4*taylor_tail(rad,cutoff-4),
                 refinement_max_absolute=err,norm=norm,curvature=curvature,normalized_curvature=curvature/norm,
                 norm_components=[ev(z) for z in norms],normalized_raw=fine[0].tolist(),
                 YU_normalized_by_sqrt_norms=fine[1].tolist())
        if ll in [.25,.1]:
            checks=[]
            for nodes in [128,192]:
                xx,wt=np.polynomial.legendre.leggauss(nodes);theta=(xx+1)*pi/2;wt=wt*pi/2
                rho=np.ones(nodes)/(2*pi);j=0
                while qq**j/(1-qq)>1e-18:
                    rho*=abs(1-qq**j*np.exp(2j*theta))**2*(1-qq**(j+1));j+=1
                energies=rad*np.cos(theta);av=(energies[:,None]+energies[None,:])/2
                om=energies[:,None]-energies[None,:];kk=bj*av**2+bd*om**2+cc*kap
                measures=[]
                for degree in [1,2,3,4]:
                    rr=qq**degree;num=1.;den=np.ones((nodes,nodes))
                    for k in range(j):
                        v=rr*qq**k;num*=1-rr*rr*qq**k
                        den*=(1-2*v*np.cos(theta[:,None]-theta[None,:])+v*v)*(1-2*v*np.cos(theta[:,None]+theta[None,:])+v*v)
                    measures.append((rho*wt)[:,None]*(rho*wt)[None,:]*num/den)
                m1,m2,m3,m4=measures
                selfmeasure=V4*m4+V2*m2+bf*bf*VF*av**2*m1+br*br*VR*av**2*m3+W*kk**2*m2
                crossmeasure=W*kk*m2
                vals=np.array([[np.sum(selfmeasure*np.cos(t*om))/norm for t in times],
                               [np.sum(crossmeasure*np.cos(t*om))/sqrt(W*norm) for t in times]])
                er=float(np.max(abs(vals-fine)));assert er<5e-10
                # Keep the ordinary-time phase cancellation visible near small values.
                cond=[float(np.sum(abs(selfmeasure*np.cos(t*om)))/max(abs(np.sum(selfmeasure*np.cos(t*om))),1e-300)) for t in times]
                checks.append(dict(nodes=nodes,normalized_values=vals.tolist(),max_absolute_error=er,
                                   mass=float(selfmeasure.sum()),phase_cancellation_condition=cond,
                                   curvature=float(np.sum(selfmeasure*om**2))))
                assert abs(np.sum(selfmeasure*om**2)/curvature-1)<5e-10
            row['independent_mehler']=checks
        records.append(row);print('raw curves',ll,'norm',norm,'passed',flush=True)
    f3=1/np.cosh(2*times)**6;f4=1/np.cosh(2*times)**8
    return dict(times=times.tolist(),records=records,small_lambda_normalized_kernel=((9*f4-2*f3)/7).tolist(),
                small_lambda_YU_normalized=(2*f3/sqrt(14)).tolist(),
                error_scope='Weighted Jacobi analytic truncation bound for raw normalized self kernel; floating arithmetic and quadrature checked by refinement, not interval arithmetic. Absolute errors used through zero crossings.')


def main():
    target=OUT/'kernel-checks.json'
    if target.exists():raise FileExistsError(target)
    a=assembly();c=curves()
    names=['code/memory_raw_kernel.py','code/memory_raw_spectrum.py','code/singlet_otoc.py',
           'results/susskind-memory-raw-spectrum/symbolic-spectrum.json',
           'results/susskind-memory-oriented-kernel/reverse-cross.json']
    target.write_text(json.dumps(dict(status='passed',assembly=a,kernel=c,
        environment=dict(python=platform.python_version(),numpy=np.__version__,sympy=s.__version__),
        input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names}),indent=2)+'\n')


if __name__=='__main__':main()
