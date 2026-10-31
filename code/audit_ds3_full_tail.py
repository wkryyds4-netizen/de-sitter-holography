"""Fraction-polynomial and Decimal-tail replay, independent of evaluator module."""
import hashlib
import json
from decimal import Decimal as D, localcontext
from fractions import Fraction as F
from math import comb
from pathlib import Path


def mul(a,b):
    out=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    return out


def power(a,n):
    out=[F(1)]
    for _ in range(n):out=mul(out,a)
    return out


def polynomial(alpha,beta):
    a=[F(1)];b=[F(1),F(2),F(1)]
    for x in alpha:a=mul(a,[F(1,4)+x*x,F(1),F(1)])
    for x in beta:b=mul(b,[1+x*x,F(2),F(1)])
    left=mul(b,power([F(1),F(1)],3));right=mul(a,power([F(2),F(1)],3))
    result=[x-y for x,y in zip(left,right)]
    while result[-1]==0:result.pop()
    return result


def audit(root):
    base=root/'results/stage4/ds3-full-tail';study=json.loads((base/'study.json').read_text())
    amp=json.loads((base/'amplification.json').read_text());errors=[];count=0
    def req(ok,msg):
        nonlocal count
        count+=1
        if not ok:errors.append(msg)
    def dc(x):return D(x.numerator)/D(x.denominator)
    def norm(z):return sum(D(x)**2 for x in z).sqrt()
    def relative(a,b,tol='1e-95'):return abs(a-b)<=D(tol)*max(abs(a),abs(b),D('1e-500'))
    refinements=[]
    with localcontext() as ctx:
        ctx.prec=180
        for row in study['rows']:
            M=F(row['M']);a=M+F(1,10);b=M+F(3,10);c=M;d=M-F(1,5);p=q=F(3,5)
            propagated=D(0)
            for sd,shadow in zip([1,-1],row['details']['shadows']):
                for sb,om in zip([1,-1],shadow['omega']):
                    alpha=[(sb*b+a+q)/2,(sb*b+c+p)/2,(sb*b+c-p)/2,(sb*b+a-q)/2]
                    hp0=(sb*b+a+c)/2;z=hp0-sd*d/2
                    groups=[(alpha,[2*hp0,sb*b,F(0)])] if z==0 else [([v-z for v in alpha],[sd*d,sb*b-z,-z]),(alpha,[2*hp0-z,sb*b,z])]
                    tail_parts=[]
                    for (aa,bb),record in zip(groups,om['details']['groups']):
                        N=record['N'];cert=record['certificate'];N0=cert['N0'];poly=polynomial(aa,bb)
                        req(len(poly)==11 and poly[-1]==1,'degree/leading coefficient')
                        req(list(reversed(poly))==list(map(F,cert['coefficients'])),'independent ratio polynomial')
                        shifted=[sum(poly[j]*comb(j,k)*N0**(j-k) for j in range(k,len(poly))) for k in range(len(poly))]
                        req(all(x>=0 for x in shifted) and list(reversed(shifted))==list(map(F,cert['shifted_coefficients'])),'exact shifted positivity')
                        req(N>=N0,'cutoff beyond certificate')
                        ad=list(map(dc,aa));bd=list(map(dc,bb));t2=D(1)
                        for n in range(N):
                            numerator=D(1);denominator=D(n+1)**2
                            for x in ad:numerator*= (D(n)+D('.5'))**2+x*x
                            for y in bd:denominator*=D(n+1)**2+y*y
                            t2*=numerator/denominator
                        t=t2.sqrt()
                        req(relative(t,norm(record['first_omitted'])),'first omitted modulus from rational product')
                        if z==0:
                            L=sum((D('.25')+(y-x)**2).sqrt() for y,x in zip(bd+[D(0)],ad))
                            req(relative(L,D(record['L'])),'digamma Lipschitz constant')
                            tail=t*L*(1+1/(D(N)+D('.5')))*(N+1)/(D(N)+D('.5'))*(1/D(N+1)+D(2)/3)
                        else:tail=t*(2*N+3)
                        tail_parts.append(tail)
                    if z==0:tail=tail_parts[0]
                    else:tail=sum(norm(co)*e for co,e in zip(om['details']['coefficients'],tail_parts))
                    req(relative(tail,D(om['tail'])),'Omega tail propagation')
                x,y=shadow['omega'];ex,ey=D(x['tail']),D(y['tail'])
                propagated+=norm(shadow['coefficient'])*(norm(x['value'])*ey+norm(y['value'])*ex+ex*ey)
            propagated*=norm(row['details']['normalization'])
            req(relative(propagated,D(row['tail'])),'full-U product and shadow tail')
            difference=[D(x)-D(y) for x,y in zip(row['value'],row['reference'])]
            req(relative(norm(difference),D(row['reference_error'])),'full complex reference residual')
            req(D(row['reference_error'])<=D(row['tail']) and row['inside_tail'],'observed residual inside bound (empirical reference)')
        for M in [2,8,32]:
            a,b=[r for r in study['rows'] if r['M']==M and r['refinement']==2]
            delta=norm([D(x)-D(y) for x,y in zip(a['value'],b['value'])])
            rel=delta/max(norm(b['value']),D('1e-500'))
            req(rel<D('1e-95'),'precision refinement agrees but does not repair truncation')
            refinements.append(dict(M=M,relative_precision_change=str(rel)))
        req(D(study['rows'][-1]['reference_error'])>D('1e43'),'retain failed high-mass approximation')
        req(amp['status']=='passed' and len(amp['rows'])==3,'elementary shadow modulus checks')
        for r in amp['rows']:
            req(max(map(D,r['relative_errors']))<D('1e-140'),'gamma vs elementary modulus')
            req(D(r['ratio_identity_error'])<D('1e-140'),'exact exponential ratio identity')
    for record in [study,amp]:
        for name,sha in record['hashes'].items():
            req(hashlib.sha256((root/name).read_bytes()).hexdigest()==sha,'input hash '+name)
    return dict(status='passed' if not errors else 'failed',checks=count,errors=errors,precision_refinements=refinements,
                scope='Exact Fraction certificates and independent Decimal replay. Does not certify Gamma/finite-sum rounding or reference values.')


if __name__=='__main__':
    root=Path(__file__).resolve().parents[1];target=root/'results/stage4/ds3-full-tail/audit.json'
    if target.exists():raise RuntimeError('Preserve prior audit')
    result=audit(root)
    result['audit_source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    target.write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2));raise SystemExit(result['status']!='passed')
