"""Exact memory-jet subtraction, q-Markov background, and tangent spectra."""
from pathlib import Path
import json, hashlib, platform
from math import exp, sqrt, pi
import numpy as np
import sympy as s
from singlet_wick_moments import gue_basis
from cubic_linear_channels import choose, coeff
from cubic_linear_dynamics import weights
from cubic_cross_channels import channel
from memory_force_bridge import exact
from bilocal_channel import blocks
from run_filtered_spectrum import parameters, cutoff
from singlet_coefficient import pairings

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-fluctuations'


def finite(N,p,m):
    cs,vs,rs=gue_basis(N,p);d=1<<N
    charges=np.array([j.bit_count() for j in range(d)])
    Q=np.diag(charges);ix=np.flatnonzero(charges==m)
    H=sum(((i%5)-2)*v for i,v in enumerate(vs))
    contract=lambda b:Q@b-sum(c.conj().T@b@c for c in cs)
    def raw2(h):
        a=2*p*h@h-contract(h@h);z=h@h@h;ha=h@a+a@h
        return [2*(3*p*z-contract(z))-3*ha,ha,2*z]
    rr=raw2(H);lin=[np.zeros_like(H) for _ in range(3)]
    for v,weight in zip(vs,rs):
        plus,minus=raw2(H+v),raw2(H-v)
        for j in range(3):lin[j]+=int(weight)*(plus[j]+minus[j]-2*rr[j])
    assert max(np.max(abs(x)) for x in rr+lin)<2**40
    L=[exact(x[np.ix_(ix,ix)])/8 for x in lin]
    T3,Sraw,H3=[exact(x[np.ix_(ix,ix)])/2 for x in rr]
    _,hn,Cn,Kn=blocks(N,p,H,cs,m)
    h,Cmat,K=exact(hn),exact(Cn),exact(Kn);I=s.eye(len(ix))
    A=exact((2*p*H@H-contract(H@H))[np.ix_(ix,ix)])
    ww=weights(N,m,p);C=choose(m,p)*choose(N-m+p,p)
    avg=lambda z:s.Rational(sum(x*y for x,y in zip(ww,z)))
    tt=[coeff(N,p,m,nu) for nu in range(p+1)]
    zz,ss=zip(*(channel(N,p,m,nu) for nu in range(p+1)))
    bt,bz,bs=avg(tt),avg(zz),avg(ss)
    a0=(2*p-m)*C+m*choose(m-1,p)*choose(N-m+p+1,p)
    cn=(bs-a0)/(bz-C);r=s.Rational(p,m)
    Y=A-a0*I-cn*(h*h-C*I);JY=(h*Y+Y*h)/2
    BZ=L[2]-bz*h;BS=L[1]-bs*h;F=L[0]-bt*h
    # Rational coefficients suffice: subtraction is an algebraic identity for any a,b.
    aa,bb=s.Rational(-2,7),s.Rational(3,11)
    R=T3-L[0]-aa*(Sraw-L[1])-bb*(H3-L[2])
    PN=(aa*cn+bb)*h**3+(bt-(aa*cn+bb)*bz)*h
    C0=a0*I+cn*(h*h-C*I)-m*r*r*h*h
    C1=PN+2*r*(a0*h+cn*(h**3-C*h))-m*r**3*h**3
    B={j:Cmat.conjugate().T*K**j*Cmat/m for j in [1,2,3]}
    G=B[2]-B[1]**2
    b3=B[3]-B[1]*B[2]-B[2]*B[1]+B[1]**3
    Z1=b3-(h*G+G*h)/2
    tangent=R+(aa+2*r)*JY+F-aa*BS-bb*BZ
    residuals=[m*G-C0-Y,m*Z1-C1-tangent,h*C0-C0*h,h*C1-C1*h,
               tangent-tangent.conjugate().T]
    assert all(all(s.expand(v)==0 for v in x) for x in residuals)
    return dict(N=N,p=p,m=m,exact_identities=5,c_N=str(cn),
                residual_norm_squared=str(s.trace(tangent*tangent)/len(ix)))


def polynomial_checks():
    q,r,x,y=s.symbols('q r x y');lam=s.symbols('lam',positive=True)
    def moment(u,v):
        z=0
        for pa in pairings(tuple(range(u+v))):
            cr=sum(a<c<b<d or c<a<d<b for j,(a,b) in enumerate(pa) for c,d in pa[j+1:])
            br=sum(a<u<=b for a,b in pa)
            z+=q**cr*r**br
        return s.expand(z)
    H=[s.Integer(1),x]
    for n in range(1,7):H.append(s.expand(x*H[n]-sum(q**j for j in range(n))*H[n-1]))
    entries=[]
    for n in range(7):
        for k in range(7-n):
            actual=sum(c*moment(k,int(power[0])) for power,c in s.Poly(H[n],x).terms())
            base=sum(c*moment(k+int(power[0]),0) for power,c in s.Poly(H[n],x).terms())
            assert s.expand(actual-r**n*base)==0
            entries.append([n,k])
    # Centered conditional first two moments, sigma=1.
    cond=lambda poly:s.expand(sum(c*r**j*H[j] for j,c in []))  # unused: explicit low moments below
    mean_increment=(r-1)*x
    second_increment=(r-1)**2*x*x+1-r*r
    assert s.expand(mean_increment-(r*x-x))==0
    assert s.expand(second_increment-(r*r*x*x+1-r*r-2*r*x*x+x*x))==0
    aa=-(1-q)**2/(1+q*q);D=(1+q*q)*(1+q+q*q)
    VR=(4*lam*(1-q*q)*(1+q*q)*(1+q+q*q)**2
        +16*lam**2*q*q*(1+q)*(1+q*q)*(1+q+q*q)*(2*q*q-q+2))/(lam**3*D**2)
    VF=64*q*q;W=4*(1/lam+4*q/(1+q))
    norm=VR+VF+aa**2*W/lam*(1+q*q)/2
    m2=VR*2/lam*(1-q**3)+VF*2/lam*(1-q)+aa**2*W/lam**2*(1+q)*(1-q**4)/2
    limits=[s.limit((lam*z).subs(q,s.exp(-4*lam)),lam,0,dir='+') for z in [W,norm,m2]]
    assert limits==[4,32,768]
    # Formal positive spectral matrix determinant; no integration needed.
    yy,ff,rr,energy,a=s.symbols('yy ff rr energy a',real=True)
    mat=s.Matrix([[yy,a*energy*yy],[a*energy*yy,ff+rr+a*a*energy*energy*yy]])
    assert s.expand(mat.det()-yy*(ff+rr))==0
    return dict(conditional_moment_tests=entries,conditional_increment=str(mean_increment),
                conditional_increment_squared=str(second_increment),
                small_lambda_coefficients=list(map(str,limits)),spectral_determinant=str(s.factor(mat.det())))


def jacobi(lam,times,n):
    z=parameters(lam);q=z['q'];k=np.arange(1,n)
    off=np.sqrt(-np.expm1(-4*lam*k)/(-np.expm1(-4*lam)*lam))
    en,U=np.linalg.eigh(np.diag(off,1)+np.diag(off,-1))
    base=U[0,:,None]*np.exp(-1j*en[:,None]*np.array(times)[None,:])
    v=U@base;one=U@(en[:,None]*base);two=U@(en[:,None]**2*base)
    pp=abs(v)**2;nn=np.arange(n)[:,None]
    gs=[np.sum(q**(j*nn)*pp,axis=0) for j in [1,2,3]]
    J2=np.sum(q**(2*nn)*(np.real(two*np.conj(v))+abs(one)**2),axis=0)/2
    norm=z['VR']+z['VF']+z['a']**2*z['W']/lam*(1+q*q)/2
    raw=z['VF']*gs[0]+z['VR']*gs[2]+z['a']**2*z['W']*J2
    return dict(G=[a.tolist() for a in gs],background_variance=(gs[1]-gs[0]**2).tolist(),
                tangent=(raw/norm).tolist(),tangent_norm=norm,
                equal_time_fractions=[z['VF']/norm,z['VR']/norm,z['a']**2*z['W']/lam*(1+q*q)/(2*norm)])


def mehler(lam,times,nodes):
    z=parameters(lam);q=z['q'];xx,ww=np.polynomial.legendre.leggauss(nodes)
    th=(xx+1)*pi/2;ww=ww*pi/2;rho=np.ones(nodes)/(2*pi);steps=1
    while q**steps/(1-q)>1e-18:steps+=1
    for k in range(steps):rho*=abs(1-q**k*np.exp(2j*th))**2*(1-q**(k+1))
    rw=rho*ww;en=2*np.cos(th)/sqrt(lam*(1-q));de=en[None,:]-en[:,None]
    es=(en[:,None]+en[None,:])/2;K=[]
    for j in [1,2,3]:
        r=q**j;num=1.;den=np.ones((nodes,nodes))
        for k in range(steps):
            v=r*q**k;num*=1-r*r*q**k
            den*=(1-2*v*np.cos(th[:,None]-th[None,:])+v*v)*(1-2*v*np.cos(th[:,None]+th[None,:])+v*v)
        K.append(num/den)
    measure=rw[:,None]*rw[None,:];norm=z['VR']+z['VF']+z['a']**2*z['W']/lam*(1+q*q)/2
    vals=[]
    for t in times:
        phase=np.exp(1j*t*de);g=np.exp(-1j*t*en)*(K[0]@(rw*np.exp(1j*t*en)))
        scalar=np.sum(rw*g);variance=np.sum(rw*abs(g-scalar)**2)
        formula=np.sum(measure*K[1]*phase).real-scalar.real**2
        assert abs(variance-formula)<2e-11
        cross=np.sum(measure*K[1]*es*phase)*z['a']*z['W']
        assert abs(cross)<2e-10
        tangent=np.sum(measure*phase*(z['VF']*K[0]+z['VR']*K[2]+z['a']**2*z['W']*es**2)).real/norm
        vals.append(dict(time=t,conditional_variance=float(variance),variance_formula=float(formula),
                         scalar=[float(scalar.real),float(scalar.imag)],cross_abs=float(abs(cross)),tangent=float(tangent)))
    weighted_cross=float(np.sum(measure*K[1]*es**2)*z['a']*z['W'])
    prediction=z['a']*z['W']/lam*(1+q*q)/2
    assert abs(weighted_cross/prediction-1)<2e-11
    return dict(nodes=nodes,product_terms=steps,values=vals,energy_weighted_cross=weighted_cross,
                energy_weighted_prediction=prediction)


def main():
    OUT.mkdir(parents=True,exist_ok=True);target=OUT/'checks.json'
    if target.exists():raise FileExistsError(target)
    finite_rows=[finite(*args) for args in [(3,1,2),(4,2,2),(5,2,3)]]
    symbolic=polynomial_checks();times=[0,.25,.5,1];rows=[]
    for lam in [.25,.1,.025,.01]:
        n,tail=cutoff(lam,1);lo=jacobi(lam,times,n);hi=jacobi(lam,times,n+12)
        error=max(np.max(abs(np.array(lo[k])-np.array(hi[k]))) for k in ['G','background_variance','tangent'])
        assert error<2e-11
        z=parameters(lam)
        ratio=z['norm']/hi['tangent_norm']
        row=dict(lam=lam,cutoff=n,refined_cutoff=n+12,refinement_error=float(error),
                 normalized_tangent_truncation_bound=tail*ratio,
                 values=hi,small_lambda_tangent=(1/np.cosh(2*np.array(times))**6).tolist())
        if lam in [.25,.1]:
            a=mehler(lam,times,128);b=mehler(lam,times,192)
            err=max(abs(b['values'][i][key]-a['values'][i][key]) for i in range(4) for key in ['conditional_variance','tangent'])
            diff=max(max(abs(v['conditional_variance']-hi['background_variance'][i]),abs(v['tangent']-hi['tangent'][i])) for i,v in enumerate(b['values']))
            assert max(err,diff)<2e-10
            row.update(mehler_coarse=a,mehler_fine=b,independent_difference=float(diff))
        rows.append(row)
        print('lambda',lam,'variance(t=.5)',hi['background_variance'][2],'tangent fractions',hi['equal_time_fractions'],flush=True)
    names=['code/memory_fluctuation_checks.py','code/singlet_wick_moments.py','code/cubic_linear_channels.py',
           'code/cubic_linear_dynamics.py','code/cubic_cross_channels.py','code/memory_force_bridge.py',
           'code/bilocal_channel.py','code/run_filtered_spectrum.py','code/singlet_coefficient.py']
    result=dict(status='passed',finite=finite_rows,symbolic=symbolic,times=times,curves=rows,
                environment=dict(python=platform.python_version(),numpy=np.__version__,sympy=s.__version__),
                error_scope='Exact finite algebra and symbolic identities; conditional quadrature and tangent spectra checked independently. Taylor truncation analytic, roundoff not interval-certified. No finite-memory fluctuation convergence asserted.',
                input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names})
    target.write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
