"""Independent saved arithmetic, elementary bound and quadrature-node replay."""
import hashlib,json
from decimal import Decimal as D,localcontext
from pathlib import Path
from fractions import Fraction as F
from math import comb
import sympy as sp
from audit_ds3_peak_limit_refined import machin


def audit(root):
    base=root/'results/stage4/ds3-uniform-local'
    study=json.loads((base/'study.json').read_text());val=json.loads((base/'validation.json').read_text());integ=json.loads((base/'integrals.json').read_text())
    errors=[];count=0
    def req(ok,msg):
        nonlocal count
        count+=1
        if not ok:errors.append(msg)
    def close(a,b,msg,tol='1e-85'):
        req(abs(a-b)<=D(tol)*max(abs(a),abs(b),D('1e-3000')),msg)
    def z(v):return tuple(D(x) for x in v)
    def norm(v):return (v[0]*v[0]+v[1]*v[1]).sqrt()
    def mul(a,b):return a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]
    def asinh(t):return (t+(1+t*t).sqrt()).ln()
    def sinh(t):return (t.exp()-(-t).exp())/2
    def cosh(t):return (t.exp()+(-t).exp())/2
    refinements=[]
    with localcontext() as ctx:
        ctx.prec=160;pi=machin()[0];gamma=D(str(sp.N(sp.EulerGamma,150)))
        req(study['status']=='passed' and len(study['rows'])==21,'21 uniform points')
        for i,row in enumerate(study['rows']):
            co=row['Ydetails']['constants'];M,u,v,R=map(D,[row['M'],row['u'],row['v'],row['R']]);B=max(abs(u),abs(v))+R;m=M-B
            C0=cosh(pi/2)*(1+R/(2*m)).sqrt()*(pi*R/4).exp()/(1-(-2*pi*m).exp()).sqrt()
            P1=lambda t:gamma+1+asinh(t)
            Q=P1(R/2)+2*P1(M+B)+2*(gamma+2*D(2).ln()+2+asinh(D(1)))
            RF=(pi*pi/6-1)/(1+m*m);RY=2*C0*((pi*pi/6-1)*Q+14)/(1+m*m)
            Smax=2*C0*Q;Ymax=Smax+RY;Dmax=C0/(4*pi)*(sinh(pi*R/2)/(pi*R/2)).sqrt()
            for name,num in dict(m=m,C0=C0,Q=Q,RF=RF,RY=RY,Smax=Smax,Ymax=Ymax,Dmax=Dmax).items():close(num,D(co[name]),f'{i}: {name}')
            N=row['Ydetails']['N'];prod=D(1)
            for j in range(1,N+1):prod*=D(j*j)/(j*j+m*m)
            ey=2*C0*prod/N*(Q+7*(2+D(N).ln()));close(ey,D(row['Ytail']),f'{i}: uniform tail')
            req(norm(z(row['D0']))<=Dmax,f'{i}: D bound')
            req(norm(z(row['Ydetails']['S']))<=Smax and norm(z(row['Y']))<=Ymax+D(row['Ytail']),f'{i}: S/Y bound')
            fin=2*mul(z(row['D0']),z(row['Ydetails']['S']))[0];close(fin,D(row['finite']),f'{i}: finite approx')
            x,y=z(row['X']),z(row['Y']);cp,cm=map(z,row['coefficients']);xy=mul(x,y)
            plus=mul(cp,xy);minus=mul(cm,(xy[0],-xy[1]));full=(plus[0]+minus[0],plus[1]+minus[1]);saved=z(row['value'])
            req(norm((full[0]-saved[0],full[1]-saved[1]))<D('1e-95'),f'{i}: complex full value')
            err=norm((saved[0]-fin,saved[1]));close(err,D(row['finite_error']),f'{i}: finite error')
            bound=2*Dmax*(RY+RF*Ymax)+2*norm(z(row['E0']))*D(row['BG'])*Ymax
            close(bound,D(row['uniform_finite_bound']),f'{i}: propagated bound');req(err<=bound+D(row['tail']),f'{i}: inequality')
            poly=list(map(F,row['Xdetails']['groups'][1]['certificate']['coefficients']));Nstar=row['Nstar'];asc=list(reversed(poly))
            req(all(sum(asc[j]*comb(j,k)*Nstar**(j-k) for j in range(k,11))>0 for k in range(11)),f'{i}: companion polynomial')
        req(val['status']=='passed' and len(val['rows'])==11,'11 independent representation checks')
        for row in val['rows']:
            old=study['rows'][row['index']];a=z(old['value']);b=z(row['reference']);diff=norm((a[0]-b[0],a[1]-b[1]))
            req(abs(diff-D(row['difference']))<D('1e-95'),f"reference difference {row['index']}")
            req(diff<=D(old['tail'])+D(row['reference_tail'])+D(row['rounding_allowance']),f"reference budget {row['index']}")
            C0,Q=map(D,[old['Ydetails']['constants'][k] for k in ['C0','Q']])
            req(max(norm(z(x)) for x in row['regularized_coefficients'])<=C0,'regularized coefficients')
            req(sum(norm(z(x)) for x in row['coefficient_derivatives'])<=C0*Q,'coefficient derivative bound')
        nodes=[json.loads(line) for line in (base/'integral-nodes.jsonl').read_text().splitlines()]
        req(len(nodes)==144 and integ['status']=='completed','144 local quadrature nodes')
        for row in integ['full_kernel']:
            subset=[x for x in nodes if x['M']==row['M'] and x['order']==row['order']]
            req(len(subset)==row['order'],'all nodes present')
            for k in range(3):
                re=sum(D(x['weight'])*D(x['test'])*D(x['r'])**k*D(x['value'][0]) for x in subset)
                im=sum(D(x['weight'])*D(x['test'])*D(x['r'])**k*D(x['value'][1]) for x in subset)
                saved=z(row['integrals'][k]);req(norm((re-saved[0],im-saved[1]))<D('1e-55'),'quadrature replay')
                tail=sum(D(x['weight'])*D(x['test'])*abs(D(x['r']))**k*D(x['tail']) for x in subset)
                close(tail,D(row['propagated_series_tails'][k]),'quadrature series tail','1e-50')
        for M in [64,256]:
            a,b=[x for x in integ['full_kernel'] if x['M']==M]
            changes=[norm(tuple(x-y for x,y in zip(z(t),z(s)))) for t,s in zip(a['integrals'],b['integrals'])]
            refinements.append(dict(M=M,absolute_changes=list(map(str,changes)),fine_series_tails=b['propagated_series_tails']))
        for record in [study,val,integ]:
            for name,sha in record['hashes'].items():req(hashlib.sha256((root/name).read_bytes()).hexdigest()==sha,name)
    x=sp.symbols('x',positive=True)
    req(sp.expand(x**3-(x*x-x+sp.Rational(1,2))*(x+1))==(x-1)/2,'base coefficient inequality')
    return dict(status='failed' if errors else 'passed',checks=count,errors=errors,quadrature_refinements=refinements,
                scope='Independent Decimal reconstruction, exact companion polynomials and all saved quadrature nodes. Gamma/Euler constants and rounding are not interval certified; numerical tests do not prove the uniform theorem.')


if __name__=='__main__':
    root=Path(__file__).resolve().parents[1];dest=root/'results/stage4/ds3-uniform-local/audit.json'
    if dest.exists():raise RuntimeError('Preserve prior audit')
    result=audit(root);result['audit_source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print(result['status'],result['checks'],result['errors'])
    raise SystemExit(result['status']!='passed')
