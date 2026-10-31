"""Independent Gaussian-return quadrature and thermal identity replay."""
import hashlib,json
from pathlib import Path
from decimal import Decimal as D,localcontext
import numpy as np
import sympy as sp
from audit_ds3_peak_limit_refined import machin

def gauss(M,p,s,t,sigma,order,limit=False):
    lo=-32.0 if limit else max(s-M,-32.0)
    edges=sorted(set([lo,32.0,0.0]+[float(j) for j in range(-32,33,2) if lo<j<32]))
    edges=[x for x in edges if x>=lo];xs,ws=np.polynomial.legendre.leggauss(order)
    u=np.concatenate([a+(b-a)*(xs+1)/2 for a,b in zip(edges,edges[1:])])
    w=np.concatenate([ws*(b-a)/2 for a,b in zip(edges,edges[1:])])
    P=np.pi;B=np.cosh(P*u);C=np.cosh(P*p)
    if limit:
        delta=u;rate=1/(8*P*(B+C))
    else:
        a=M+u;sm=np.exp(-P*(a+M));invA=2*sm/(1+sm*sm)
        bal=(1-(B+C)*invA/(1+C*invA))/(8*P*(B+C))
        g=lambda x:x*np.exp(P/(np.sqrt(1+x*x)+x))/(1+np.exp(-2*P*x))
        delta=u*(a+M)/(np.sqrt(1+a*a)+np.sqrt(1+M*M));rate=g(a)/g(M)*bal
    return complex(np.sum(w*rate*np.exp(-sigma*sigma*delta**2-P*delta-1j*t*delta)))

def audit(root):
    base=root/'results/stage4/ds3-fock-vacuum';data=json.loads((base/'gaussian.json').read_text())
    count=0;errors=[];replay=[]
    def req(ok,name):
        nonlocal count
        count+=1
        if not ok:errors.append(name)
    req(data['status']=='passed' and len(data['rows'])==16 and len(data['limits'])==4,'Gaussian inventory')
    with localcontext() as ctx:
        ctx.prec=100;pi=machin()[0];sig=D(data['sigma']);s=D('.2')
        Lam=(pi*pi/(4*sig*sig)).exp();Bs=2*(pi/((1+s*s).sqrt()+s)).exp()
        req(abs(Lam-D(data['multiplier_bound']))<D('1e-54')*Lam,'multiplier bound')
        for row in data['rows']:
            x=row['fine'];M,p,t=map(float,[x[k] for k in ['M','p','t']]);sigma=float(sig)
            q1=gauss(M,p,.2,t,sigma,96);q2=gauss(M,p,.2,t,sigma,160);ref=complex(*map(float,x['value']))
            req(abs(q2-ref)<1e-13 and abs(q2-q1)<1e-13,'independent Gaussian quadrature')
            re,im=map(D,x['value'])
            req(abs(4*im*im-D(x['vacuum_commutator_norm']))<D('1e-56'),'vacuum norm')
            for y in [row['coarse'],row['fine']]:
                cut=D(y['D']);tail=Lam*Bs*(-pi*cut).exp()/(2*pi)*(1/pi+(cut/pi+1/pi**2)/s)
                req(abs(tail-D(y['tail']))<D('1e-53')*tail,'Gaussian spectral tail')
            ca,cb=map(D,row['coarse']['value']);diff=((re-ca)**2+(im-cb)**2).sqrt()
            req(abs(diff-D(row['difference']))<D('1e-56'),'cutoff difference')
            req(diff<=D(row['coarse']['tail'])+D(x['tail'])+D(row['rounding_allowance']),'cutoff criterion')
            replay.append(dict(M=M,t=t,value=[q2.real,q2.imag],difference=abs(q2-ref),refinement=abs(q2-q1),tolerance=1e-13))
        for row in data['limits']:
            t=float(row['t']);q=gauss(0,float(data['p']),.2,t,float(sig),160,True)
            ref=complex(*map(float,row['value']))
            req(abs(q-ref)<1e-13,'independent limiting transform')
            continued=list(map(D,row['kms_continued']));reverse=list(map(D,row['reversed']))
            diff=sum((a-b)**2 for a,b in zip(continued,reverse)).sqrt()
            req(diff<D('1e-50'),'saved KMS identity')
        for row in data['detailed_balance']:req(D(row['error'])<D('1e-55'),'saved detailed balance')
        for name,sha in data['hashes'].items():req(hashlib.sha256((root/name).read_bytes()).hexdigest()==sha,name)
    x,sigma=sp.symbols('x sigma',real=True)
    req(sp.expand(-sigma**2*x*x-sp.pi*x + (sigma*x+sp.pi/(2*sigma))**2-sp.pi**2/(4*sigma**2))==0,'completing square')
    return dict(status='passed' if not errors else 'failed',checks=count,errors=errors,independent_quadratures=replay,
        scope='Independent Gaussian spectral arithmetic and quadrature; no interval rounding or full fourpoint computation.')

if __name__=='__main__':
    root=Path(__file__).resolve().parents[1];dest=root/'results/stage4/ds3-fock-vacuum/gaussian-audit.json'
    if dest.exists():raise RuntimeError('Preserve audit')
    d=audit(root);d['audit_source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    dest.write_text(json.dumps(d,indent=2),encoding='utf-8');print(d['status'],d['checks'],d['errors'])
    if d['status']!='passed':raise SystemExit(1)
