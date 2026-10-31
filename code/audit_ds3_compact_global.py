"""Independent Decimal/Fraction replay; no import of the evaluated kernel."""
import hashlib,json
from decimal import Decimal as D,localcontext
from fractions import Fraction as F
from math import comb
from pathlib import Path
import sympy as sp
from audit_ds3_peak_limit_refined import machin

def audit(root):
    base=root/'results/stage4/ds3-compact-global'
    study=json.loads((base/'study.json').read_text());repair=json.loads((base/'collision-diagnosis.json').read_text())
    count=0;errors=[]
    def req(ok,msg):
        nonlocal count
        count+=1
        if not ok:errors.append(msg)
    def close(a,b,msg,tol='1e-78'):
        req(abs(a-b)<=D(tol)*max(abs(a),abs(b),D('1e-3000')),msg)
    def z(v):return tuple(map(D,v))
    def norm(v):return sum(x*x for x in v).sqrt()
    def mul(a,b):return a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]
    def asinh(t):return (t+(1+t*t).sqrt()).ln()
    def sinh(t):return (t.exp()-(-t).exp())/2
    def cosh(t):return (t.exp()+(-t).exp())/2
    with localcontext() as ctx:
        ctx.prec=160;pi=machin()[0];gamma=D(str(sp.N(sp.EulerGamma,150)))
        req(study['status']=='failed' and len(study['rows'])==18,'retain initial failure and 18 rows')
        req([r['index'] for r in study['references'] if not r['passed']]==[16],'exactly reference 16 failed')
        for i,row in enumerate(study['rows']):
            co=row['Ydetails']['constants'];M,u,v,p,q,R=map(D,[row[k] for k in ['M','u','v','p','q','R']])
            B=max(abs(u),abs(v))+R;m=M-B;astar=max(abs(x) for x in [(v+q)/2,(p-u)/2,(-u-p)/2,(v-q)/2])+R/2
            Cseries=cosh(pi*astar)**2
            C0=cosh(pi*astar)*(1+R/(2*m)).sqrt()*(pi*R/4).exp()/(1-(-2*pi*m).exp()).sqrt()
            P1=lambda t:gamma+1+asinh(t)
            Q=P1(R/2)+2*P1(M+B)+2*(gamma+2*D(2).ln()+2+asinh(2*astar))
            RF=Cseries*(pi*pi/6-1)/(1+m*m);RY=2*C0*Cseries*((pi*pi/6-1)*Q+14)/(1+m*m)
            Smax=2*C0*Q;Ymax=Smax+RY;Dmax=C0/(4*pi)*(sinh(pi*R/2)/(pi*R/2)).sqrt()
            for name,num in dict(m=m,astar=astar,Cseries=Cseries,C0=C0,Q=Q,RF=RF,RY=RY,Smax=Smax,Ymax=Ymax,Dmax=Dmax).items():close(num,D(co[name]),f'{i}: {name}')
            N=row['Ydetails']['N'];prod=D(1)
            for j in range(1,N+1):prod*=D(j*j)/(j*j+m*m)
            ey=2*C0*Cseries*prod/N*(Q+7*(2+D(N).ln()));close(ey,D(row['Ytail']),f'{i}: tail')
            req(norm(z(row['D0']))<=Dmax,f'{i}: D bound')
            req(norm(z(row['Ydetails']['S']))<=Smax and norm(z(row['Y']))<=Ymax+D(row['Ytail']),f'{i}: S/Y bound')
            fin=2*mul(z(row['D0']),z(row['Ydetails']['S']))[0];close(fin,D(row['finite']),f'{i}: finite')
            xy=mul(z(row['X']),z(row['Y']));cp,cm=map(z,row['coefficients'])
            plus=mul(cp,xy);minus=mul(cm,(xy[0],-xy[1]));saved=z(row['value'])
            req(norm((plus[0]+minus[0]-saved[0],plus[1]+minus[1]-saved[1]))<D('1e-85'),f'{i}: full complex value')
            err=norm((saved[0]-fin,saved[1]));close(err,D(row['finite_error']),f'{i}: error')
            bound=2*Dmax*(RY+RF*Ymax)+2*norm(z(row['E0']))*D(row['BG'])*Ymax
            close(bound,D(row['uniform_finite_bound']),f'{i}: bound');req(err<=bound+D(row['tail']),f'{i}: inequality')
            poly=list(map(F,row['Xdetails']['groups'][1]['certificate']['coefficients']));asc=list(reversed(poly));Ns=row['Nstar']
            req(all(sum(asc[j]*comb(j,k)*Ns**(j-k) for j in range(k,11))>0 for k in range(11)),f'{i}: polynomial derivatives')
        for row in study['references']:
            old=z(study['rows'][row['index']]['value']);ref=z(row['value']);diff=norm(tuple(x-y for x,y in zip(old,ref)))
            req(abs(diff-D(row['difference']))<D('1e-100'),'reference difference replay')
            req((diff<=D(row['tail'])+D(study['rows'][row['index']]['tail'])+D(row['rounding_allowance']))==row['passed'],'reference criterion replay')
        for row in repair['rows']:
            old=z(study['rows'][16]['value']);ref=z(row['value']);diff=norm(tuple(x-y for x,y in zip(old,ref)))
            req(abs(diff-D(row['difference']))<D('1e-110'),'repair difference replay')
            req((diff<=D(row['tail'])+D(study['rows'][16]['tail'])+D(repair['rounding_allowance']))==row['passed'],'repair criterion')
        req(not repair['rows'][0]['passed'] and all(r['passed'] for r in repair['rows'][1:]),'failed reproduction then three successful controls')
        rr=list(map(F,repair['rounded_labels']));ee=list(map(F,repair['exact_labels']))
        req(rr[2]+rr[4]-rr[0]-rr[1]==F(repair['exact_collision_residual'])!=0,'nonzero exact rounded residual')
        req(ee[2]+ee[4]-ee[0]-ee[1]==0,'exact collision retained')
        for row in study['tails']:
            K,s=map(D,[row['K'],row['s']]);B=2*(pi/((1+s*s).sqrt()+s)).exp()
            eps=B*(-pi*K).exp()/(4*pi*pi)*(K/pi+1/pi**2+(K*K/pi+2*K/pi**2+2/pi**3)/s)
            nb=B/(4*pi**4)+B/(2*s*pi**5)
            for name,num in [('B',B),('epsilon',eps),('norm_bound',nb),('pairing_tail',2*(eps*nb).sqrt())]:close(num,D(row[name]),f'tail K{K} {name}')
        for row in study['product_checks']:
            a=D(row['alpha']);n=row['n'];sq=D(1)
            for k in range(n):sq*=((D(k)+D('.5'))**2+a*a)/D((k+1)**2)
            close(sq.sqrt(),D(row['actual']),'Pochhammer modulus');req(sq<=cosh(pi*a)/(n+1),'Pochhammer bound')
        for record in [study,repair]:
            for name,sha in record['hashes'].items():req(hashlib.sha256((root/name).read_bytes()).hexdigest()==sha,name)
    n,x,s,P=sp.symbols('n x s P',positive=True)
    req(sp.factor((n+1)**3-(n+sp.Rational(1,2))**2*(n+2))==(3*n+2)/4,'central coefficient induction')
    tail=sp.exp(-P*x)*(x/P+1/P**2+(x*x/P+2*x/P**2+2/P**3)/s)
    req(sp.simplify(sp.diff(tail,x)+x*(1+x/s)*sp.exp(-P*x))==0,'radial tail antiderivative')
    return dict(status='passed' if not errors else 'failed',checks=count,errors=errors,
        scope='Independent saved arithmetic, exact rational collision and polynomial checks, tail antiderivative. Not a proof checker, interval certification or full physical quadrature.')

if __name__=='__main__':
    root=Path(__file__).resolve().parents[1];dest=root/'results/stage4/ds3-compact-global/audit.json'
    if dest.exists():raise RuntimeError('Preserve audit')
    result=audit(root);result['audit_source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print(result)
    raise SystemExit(result['status']!='passed')
