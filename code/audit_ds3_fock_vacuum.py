"""Independent elementary quadrature, saved arithmetic and phase/path audit."""
import hashlib,json
from decimal import Decimal as D,localcontext
from pathlib import Path
import numpy as np
import sympy as sp
from audit_ds3_peak_limit_refined import machin

def quadrature(M,p,s,t,order):
    # Independent vectorized Gauss calculation; no production kernel imports.
    lo=max(s-M,-32.0);edges=sorted(set([lo,32.0,0.0]+[float(j) for j in range(-32,33,2) if lo<j<32]))
    edges=[x for x in edges if x>=lo];xs,ws=np.polynomial.legendre.leggauss(order)
    us=np.concatenate([a+(b-a)*(xs+1)/2 for a,b in zip(edges,edges[1:])])
    weights=np.concatenate([ws*(b-a)/2 for a,b in zip(edges,edges[1:])]);a=M+us;P=np.pi
    B=np.cosh(P*us);C=np.cosh(P*p);small=np.exp(-P*(a+M));invA=2*small/(1+small**2)
    bal=(1-(B+C)*invA/(1+C*invA))/(8*P*(B+C))
    ga=a*np.exp(P/(np.sqrt(1+a*a)+a))/(1+np.exp(-2*P*a))
    gb=M*np.exp(P/(np.sqrt(1+M*M)+M))/(1+np.exp(-2*P*M))
    delta=us*(a+M)/(np.sqrt(1+a*a)+np.sqrt(1+M*M))
    return complex(np.sum(weights*ga/gb*bal*np.exp(-1j*t*delta)))

def audit(root):
    base=root/'results/stage4/ds3-fock-vacuum';data=json.loads((base/'diagnostics-recovered.json').read_text())
    count=0;errors=[];replays=[];asymptotics=[]
    def req(ok,name):
        nonlocal count
        count+=1
        if not ok:errors.append(name)
    req(data['status']=='passed' and len(data['vacuum'])==20,'20 quadrature pairs')
    with localcontext() as ctx:
        ctx.prec=100;pi=machin()[0];s=D('.2');Bs=2*(pi/((1+s*s).sqrt()+s)).exp()
        bound=Bs*(1/(2*pi*pi)+1/(2*s*pi**3))
        req(abs(bound-D(data['block_bound']))<D('1e-55'),'uniform block constant')
        for row in data['rates']:
            a,b,p=map(D,[row[k] for k in ['a','b','p']])
            cosh=lambda x:(x.exp()+(-x).exp())/2
            ga=lambda x:x*(pi/((1+x*x).sqrt()+x)).exp()/(1+(-2*pi*x).exp())
            A=cosh(pi*(a+b));B=cosh(pi*(a-b));C=cosh(pi*p)
            val=ga(a)/ga(b)*(A-B)/(8*pi*(A+C)*(B+C))
            req(abs(val-D(row['rate']))<D('1e-52')*val,'independent rate')
            req(abs(D(row['raw'])-val)<D('1e-48')*val,'raw E/w identity')
        for i,row in enumerate(data['vacuum']):
            fine=row['fine'];coarse=row['coarse'];re,im=map(D,fine['value']);M,t,p=map(float,[fine[k] for k in ['M','t','p']])
            req(abs(4*im*im-D(fine['vacuum_commutator_norm']))<D('1e-58'),'vacuum norm replay')
            for record in [coarse,fine]:
                cutoff=D(record['D']);tail=Bs*(-pi*cutoff).exp()/(2*pi)*(1/pi+(cutoff/pi+1/pi**2)/s)
                req(abs(tail-D(record['tail']))<D('1e-54')*tail,'one-step tail')
            a,b=map(D,coarse['value']);diff=((re-a)**2+(im-b)**2).sqrt()
            req(abs(diff-D(row['difference']))<D('1e-58'),'cutoff difference')
            req(diff<=D(coarse['tail'])+D(fine['tail'])+D(row['rounding_allowance']),'tail/rounding criterion')
            q1=quadrature(M,p,.2,t,96);q2=quadrature(M,p,.2,t,160);reference=complex(float(re),float(im))
            req(abs(q2-reference)<1e-14 and abs(q1-q2)<1e-14,'independent Gauss comparison')
            replays.append(dict(index=i,coarse_order=96,fine_order=160,value=[q2.real,q2.imag],difference=abs(q2-reference),refinement=abs(q1-q2),tolerance=1e-14))
            first=complex(*map(float,fine['first_correction']));derivative=first.imag*M
            second=-(t+1j*np.pi)*derivative/2
            scaled=complex(*map(float,fine['second_scaled']))
            asymptotics.append(dict(M=M,t=t,p=p,second_coefficient=[second.real,second.imag],
                scaled_residual=[scaled.real,scaled.imag],third_scaled_error=[M*(scaled-second).real,M*(scaled-second).imag]))
        for r in data['finite_fock']['fourpoint']:
            a,b,c=[complex(*r[k]) for k in ['full','vacuum','two']]
            req(abs(a-b-c)<1e-14,'finite fourpoint paths')
        for r in data['finite_fock']['commutators']:req(abs(r['full']-r['vacuum']-r['two'])<1e-14,'finite orthogonal sectors')
        for r in data['finite_fock']['bounds']:req(r['norm']<=r['bound']+1e-14,'finite sqrt N bound')
        for name,sha in data['hashes'].items():req(hashlib.sha256((root/name).read_bytes()).hexdigest()==sha,name)
    u,v,r,M=sp.symbols('u v r M',real=True);ts=sp.symbols('t1:5',real=True);t1,t2,t3,t4=ts
    bra=v*t2+u*t1;ket_s=v*t3+u*t4;ket_t=(u-r)*t3+(v+r)*t4
    req(sp.expand(bra-ket_s+v*(t3-t2)+u*(t4-t1))==0,'direct phase')
    req(sp.expand((bra-ket_t).subs(r,0)+u*(t3-t1)+v*(t4-t2))==0,'exchange phase')
    e,x,t=sp.symbols('e x t',real=True)
    ratio=1+x*e-sp.pi*x*e**2/2
    phase=1+sp.I*t*x*e**2/2
    req(sp.expand(ratio*phase).coeff(e,2)==x*(-sp.pi+sp.I*t)/2,'second-order amplitude phase')
    paths=[]
    for signs in [(a,b,c,d) for a in [1,-1] for b in [1,-1] for c in [1,-1] for d in [1,-1]]:
        seq=np.cumsum(signs)
        if min(seq)>=0 and seq[-1]==0:paths.append(signs)
    req(set(paths)=={(1,1,-1,-1),(1,-1,1,-1)},'exact vacuum particle paths')
    return dict(status='passed' if not errors else 'failed',checks=count,errors=errors,independent_quadratures=replays,
        second_order_observations=asymptotics,
        scope='Independent Decimal and Gauss checks, exact phase/path identities. No interval error, full spin enumeration or all-N proof verification.')

if __name__=='__main__':
    root=Path(__file__).resolve().parents[1];dest=root/'results/stage4/ds3-fock-vacuum/audit.json'
    if dest.exists():raise RuntimeError('Preserve audit')
    data=audit(root);data['audit_source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    dest.write_text(json.dumps(data,indent=2),encoding='utf-8');print(data['status'],data['checks'],data['errors'])
    if data['status']!='passed':raise SystemExit(1)
