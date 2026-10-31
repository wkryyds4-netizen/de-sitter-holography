"""Saved-arithmetic and provenance audit, including explicit failed numerical settings."""
import hashlib
import json
from decimal import Decimal as D,localcontext
from audit_ds3_crossing_dictionary import PI


def audit(root):
    base=root/'results/stage4/ds3-collision'
    read=lambda n:json.loads((base/n).read_text(encoding='utf-8'))
    names=['study.json','validation.json','asymptotic.json','controlled.json','controlled-validation.json']
    records={n:read(n) for n in names};errors=[]
    def req(ok,msg):
        if not ok:errors.append('collision: '+msg)
    def z(v):return tuple(D(x) for x in v)
    def norm(v):return (v[0]*v[0]+v[1]*v[1]).sqrt()
    def err(v,w):return norm((v[0]-w[0],v[1]-w[1]))
    def close(a,b,msg,tol='1e-50'):req(abs(a-b)<D(tol)*max(D(1),abs(b)),msg)
    def asinh(x):return (x+(x*x+1).sqrt()).ln()
    with localcontext() as ctx:
        ctx.prec=100
        study=records['study.json'];v=records['validation.json'];asym=records['asymptotic.json'];co=records['controlled.json'];cv=records['controlled-validation.json']
        req(study['status']=='completed' and len(study['rows'])==12,'12 preserved prototype study points')
        for r in study['rows']:
            a,b=[z(k['value']) for k in r['runs']]
            close(err(a,b),D(r['absolute_difference']),'prototype precision difference','1e-45')
            close(err(a,b)/norm(b),D(r['relative_difference']),'prototype relative difference','1e-45')
        req(v['status']=='completed' and len(v['checks'])==9,'nine cross-representation checks')
        for r in v['checks']:
            e=err(z(r['value']),z(r['reference']));rel=e/max(norm(z(r['reference'])),D('1e-60'))
            close(e,D(r['absolute_error']),'check absolute','1e-65')
            close(rel,D(r['relative_error']),'check relative','1e-64')
            req(r['passed'] and rel<D(r['tolerance']),'check pass '+r['name'])
        pb=v['prototype_bias']
        e=err(z(pb['old']),z(pb['reference']));e2=err(z(pb['corrected']),z(pb['reference']))
        close(e,D(pb['old_error']),'old coefficient bias','1e-65');close(e2,D(pb['corrected_error']),'corrected coefficient error','1e-65')
        req(e>D('1e-28') and e2<D('1e-50'),'preserved shared bias')
        req(len(v['heavy'])==3,'heavy benchmark groups')
        for r in v['heavy']:
            reference=z(r['reference'])
            for a in r['attempts']:
                relative=err(z(a['value']),reference)/norm(reference)
                close(relative,D(a['relative_error']),'heavy attempted errors','1e-55')
        h32=v['heavy'][-1]
        req(h32['M']==32 and len(h32['attempts'])==5,'all M32 attempts')
        req(D(h32['attempts'][0]['relative_error'])>D('1e24'),'retained divergent short tail')
        req(D(h32['attempts'][1]['relative_error'])>D('1e6'),'retained inadequate refinement')
        req(D(h32['attempts'][-1]['relative_error'])<D('1e-23'),'selected final M32 check')
        req(asym['status']=='passed' and len(asym['rows'])==12,'analytic-window samples')
        for r in asym['rows']:
            b,d=D(r['details']['b']),D(r['details']['d'])
            tops=[z(x) for x in r['details']['top']]
            B=2+asinh(b)+asinh(d)+4*D(2).ln()*sum(norm((a-1,c)) for a,c in tops)
            bound=(PI*PI/6-1)*B/(1+min(b,d)**2)
            close(bound,D(r['bound']),'explicit leading bound','1e-60')
            e=err(z(r['value']),z(r['leading']))
            close(e,D(r['absolute_error']),'leading error','1e-60')
            req(e<bound and r['passed'],'window bound check')
        req(co['status']=='completed' and len(co['rows'])==10,'controlled direct sums')
        for r in co['rows']:
            M=D(r['M']);sc=z(r['scaled']);value=z(r['value'])
            close(value[0]/(2*M.ln()),sc[0],'scaled real','1e-85')
            close(value[1]/(2*M.ln()),sc[1],'scaled imag','1e-85')
            target=D(r['sinc_target'])
            close(norm((sc[0]-target,sc[1])),D(r['scaled_error']),'slow scaling error','1e-85')
            if r['lambda_value']==0:
                b=M+D('.3');d=M-D('.2');m=min(b,d)
                B=2+asinh(b)+asinh(d)+4*D(2).ln()*sum((D('.25')+a*a).sqrt() for a in map(D,['.2','.15','-.45','-.4']))
                product=D(1)
                for j in range(1,129):product*=D(j*j)/(D(j*j)+m*m)
                bound=B*product/128
                saved=D(r['details'][-1]['series_tail_bound'])
                req(abs(bound-saved)<=abs(bound)*D('1e-85'),'direct-sum exact-arithmetic tail')
        req(cv['status']=='passed' and len(cv['rows'])==12,'independent direct-sum comparisons')
        for r in cv['rows']:
            e=err(z(r['value']),z(r['reference']))
            close(e,D(r['error']),'independent controlled error','1e-90')
            req(r['tail_passed'] and e<=D(r['tail'])+D('1e-85'),'series tail plus declared rounding allowance')
            req(r['leading_passed'] and D(r['leading_error'])<=D(r['leading_bound']),'offsurface leading bound')
    for record in records.values():
        for path,sha in record['hashes'].items():
            req(hashlib.sha256((root/path).read_bytes()).hexdigest()==sha,'input hash '+path)
    return errors


if __name__=='__main__':
    from pathlib import Path
    errors=audit(Path(__file__).resolve().parents[1]);print(json.dumps(dict(status='failed' if errors else 'passed',errors=errors),indent=2))
    raise SystemExit(bool(errors))
