"""Independent saved-value arithmetic audit; retained failures are mandatory."""
import hashlib,json
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal as D,localcontext
from math import comb
from audit_ds3_full_tail import polynomial


def audit(root):
    base=root/'results/stage4/ds3-shadow-duality'
    read=lambda name:json.loads((base/name).read_text())
    study,validation,refinement=map(read,['study.json','validation.json','refinement.json'])
    errors=[];count=0;precision=[]
    def req(ok,msg):
        nonlocal count
        count+=1
        if not ok:errors.append(msg)
    def norm(z):return sum(D(v)**2 for v in z).sqrt()
    def diff(a,b):return norm([D(x)-D(y) for x,y in zip(a,b)])
    def dc(x):return D(x.numerator)/D(x.denominator)
    def close(a,b):return abs(a-b)<=D('1e-60')*max(abs(a),abs(b),D('1e-500'))
    with localcontext() as ctx:
        ctx.prec=170
        req(len(study['rows'])==18,'all 18 runs')
        for row in study['rows']:
            M=F(row['M']);a=M+F(1,10);b=M+F(3,10);c=M;d=M-F(1,5);p=q=F(3,5)
            for sb,om in zip([1,-1],row['details']['omega']):
                A=[(sb*b+a+q)/2,(sb*b+c+p)/2,(sb*b+c-p)/2,(sb*b+a-q)/2];hp0=(sb*b+a+c)/2;z=hp0-d/2
                groups=[(A,[2*hp0,sb*b,F(0)])] if z==0 else [([v-z for v in A],[d,sb*b-z,-z]),(A,[2*hp0-z,sb*b,z])]
                tails=[]
                for (aa,bb),record in zip(groups,om['details']['groups']):
                    cert=record['certificate'];N=record['N'];N0=cert['N0'];poly=polynomial(aa,bb)
                    shifted=[sum(poly[j]*comb(j,k)*N0**(j-k) for j in range(k,len(poly))) for k in range(len(poly))]
                    req(list(reversed(poly))==list(map(F,cert['coefficients'])),'ratio polynomial')
                    req(N>=N0 and all(x>=0 for x in shifted) and list(reversed(shifted))==list(map(F,cert['shifted_coefficients'])),'exact tail certificate')
                    ad,bd=list(map(dc,aa)),list(map(dc,bb));t2=D(1)
                    for n in range(N):
                        nu=D(1);de=D(n+1)**2
                        for x in ad:nu*=(D(n)+D('.5'))**2+x*x
                        for y in bd:de*=D(n+1)**2+y*y
                        t2*=nu/de
                    t=t2.sqrt();req(close(t,norm(record['first_omitted'])),'omitted term modulus')
                    if z==0:
                        L=sum((D('.25')+(y-x)**2).sqrt() for y,x in zip(bd+[D(0)],ad))
                        req(close(L,D(record['L'])),'digamma bound constant')
                        tail=t*L*(1+1/(D(N)+D('.5')))*(N+1)/(D(N)+D('.5'))*(1/D(N+1)+D(2)/3)
                    else:tail=t*(2*N+3)
                    tails.append(tail)
                tail=tails[0] if z==0 else sum(norm(co)*e for co,e in zip(om['details']['coefficients'],tails))
                req(close(tail,D(om['tail'])),'Omega tail')
            co=row['details']['coefficients'];x,y=row['details']['omega'];ex,ey=D(x['tail']),D(y['tail'])
            bound=(norm(co[0])+norm(co[1]))*(norm(x['value'])*ey+norm(y['value'])*ex+ex*ey)
            req(close(bound,D(row['tail'])),'reflected full tail')
            req(abs(norm(co[0])/norm(co[1])-1)<D('1e-60'),'balanced coefficient moduli')
            if row['reference'] is not None:
                e=diff(row['value'],row['reference'])
                req(abs(e-D(row['reference_error']))<D('1e-60'),'complete-kernel complex residual (saved precision floor)')
                req(e<D(row['tail']),'empirical old reference within analytic tail')
        failed=[r for r in validation['checks'] if not r['passed']]
        req(validation['status']=='failed' and len(failed)==2,'retain two original overstrict failures')
        req(len(validation['checks'])==33 and all(validation['symbolic']['checks'].values()),'original count and exact symbolic dictionary')
        for r in validation['checks']:
            e=diff(r['value'],r['reference'])/norm(r['reference'])
            allowance=D('1e-30') if r['name']=='independent_Euler_pair' else D('1e-45') if not r['name'].startswith('saved_') else D('1e-75')
            req(abs(e-D(r['relative_error']))<allowance,'original reflection relative residual at saved precision')
            req((e<D(r['tolerance']))==r['passed'],'original pass/fail classification')
        req(refinement['status']=='passed' and len(refinement['rows'])==4,'four fresh original-formula refinements')
        for r in refinement['rows']:
            e=diff(r['direct'],r['reflected'])/norm(r['direct'])
            req(abs(e-D(r['relative_error']))<D('1e-90') and e<D(r['tolerance']) and r['passed'],'fresh reflection check')
        for M in [2,8,16,32,64,128]:
            a,b=[r for r in study['rows'] if r['M']==M and r['refinement']==2]
            delta=diff(a['value'],b['value']);relative=delta/norm(b['value'])
            req(relative<D('1e-60'),'70/110-digit relative refinement')
            precision.append(dict(M=M,absolute_change=str(delta),relative_change=str(relative),
                                  final_analytic_tail=b['tail'],rounding_exceeds_tail=delta>D(b['tail'])))
        req(precision[-1]['rounding_exceeds_tail'],'retain arithmetic floor at M128')
    for data in [study,validation,refinement]:
        for path,sha in data['hashes'].items():req(hashlib.sha256((root/path).read_bytes()).hexdigest()==sha,'hash '+path)
    return dict(status='passed' if not errors else 'failed',checks=count,errors=errors,precision_refinements=precision,
                retained_original_failures=2,scope='Independent polynomial/tail/residual replay; no interval rounding, uniform mass limit, or original-reference certification.')


if __name__=='__main__':
    root=Path(__file__).resolve().parents[1];dest=root/'results/stage4/ds3-shadow-duality/audit.json'
    if dest.exists():raise RuntimeError('Preserve previous audit')
    data=audit(root);data['audit_source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    dest.write_text(json.dumps(data,indent=2),encoding='utf-8');print(json.dumps(data,indent=2));raise SystemExit(data['status']!='passed')
