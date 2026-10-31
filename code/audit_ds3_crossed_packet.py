"""Independent saved-arithmetic and provenance audit; no expensive 6j rerun."""
import hashlib
import json
from decimal import Decimal as D,localcontext
from fractions import Fraction
from math import comb
from audit_ds3_crossing_dictionary import PI


def audit(root,number,norm):
    base=root/'results/stage4/ds3-crossed-packet'
    read=lambda n:json.loads((base/n).read_text(encoding='utf-8'))
    summary=read('summary.json');bound=read('bounds-validation.json');accelerator=read('accelerator-validation.json')
    errors=[]
    def req(ok,msg):
        if not ok:errors.append('crossed packet: '+msg)
    def close(a,b,msg,tol='1e-28'):
        req(abs(a-b)<=D(tol)*max(abs(b),D(1)),msg)
    def f(x):return D.from_float(x) if isinstance(x,float) else D(x)
    def z(x):return tuple(map(f,x))
    def sub(x,y):return x[0]-y[0],x[1]-y[1]
    def ch(x):return (x.exp()+(-x).exp())/2
    def sh(x):return (x.exp()-(-x).exp())/2
    req(Fraction(630)*sum(Fraction((-1)**j*comb(4,j),j+5) for j in range(5))==1,'exact Beta probability')
    req(summary['status']=='passed' and len(summary['rows'])==5,'five completed raw runs')
    req(bound['status']=='passed' and len(bound['checks'])==24,'bounds check count')
    req(accelerator['status']=='completed' and len(accelerator['rows'])==13,'precision comparison count')
    allrecords=[summary,bound,accelerator]
    for row in summary['rows']:
        raw=read(row['label']+'.json');allrecords.append(raw)
        req(raw['status']=='completed','completed '+row['label'])
        req(hashlib.sha256((root/raw['raw_grid']).read_bytes()).hexdigest()==raw['raw_grid_sha256'],'raw grid '+row['label'])
        req(raw['points']==row['points'],'raw point count')
    for record in allrecords:
        for path,sha in record['hashes'].items():
            req(hashlib.sha256((root/path).read_bytes()).hexdigest()==sha,'input hash '+path)
    with localcontext() as ctx:
        ctx.prec=75
        for check in bound['checks']:
            obs,ex=number(check['observed']),number(check['expected'])
            error=norm(sub(obs,ex))/max(norm(ex),D('1e-30'))
            req(check['passed'] is True and error<D(check['tolerance']),'saved '+check['name'])
            close(error,D(check['error']),'saved relative error','1e-30')
        s,c0,c1,p,q=map(D,['.2','.4','.8','.6','.9'])
        delta=lambda x:1/((1+x*x).sqrt()+x)
        Cs=(PI*delta(s)).exp();C0=(Cs-1)/(2*delta(s))
        gmin=c0*(PI*delta(c1)).exp()/(1+(-2*PI*c0).exp())
        ap=p/(4*PI*sh(PI*p));aq=q/(4*PI*sh(PI*q))
        nmax=((c1+C0)*ap*aq+(ap+aq)/(2*PI**3))/gmin
        for row in bound['tails']:
            L=D(row['L']);decay=(-PI*(L-c1)).exp()
            ta=(PI*delta(L)).exp()*decay/(16*PI**3*gmin)*(L*L+L*(3/PI-c0)+3/PI**2-c0/PI)
            tb=decay/(8*PI**4*gmin)*(L+C0+2/PI)
            for key,value in [('gmin',gmin),('C0',C0),('tail_a',ta),('tail_b',tb),
                              ('tail_norm_bound',ta+tb),('norm_squared_bound',nmax),
                              ('crossed_absolute_bound',2*((ta+tb)*nmax).sqrt())]:
                close(value,D(row[key]),'analytic '+key)
        ns,nt=map(f,[summary['Ns'],summary['Nt']])
        for i,value in enumerate(summary['X']):
            re,im=z(value);diff=ns+nt-2*re
            close(diff,f(summary['D'][i]),'order difference','1e-19')
            ratio=norm((re,im))/(ns*nt).sqrt()
            close(ratio,f(summary['cauchy_ratios'][i]),'Cauchy ratio','1e-14')
            req(ratio<=1 and diff>0,'positive physical overlap')
        w=bound['witnesses'][-1]
        witness=f(w['overlap'])**2/f(w['test_norm_squared'])
        close(witness,f(w['order_difference_lower_bound']),'witness Cauchy bound','1e-23')
        req(f(w['overlap'])<0 and f(w['test_norm_squared'])>0 and
            ch(PI*p)*ch(PI*q)>sh(PI*s)**2,'strict witness sign')
        req(f(summary['D'][0])>witness,'full overlap consistent with witness')
        bad=[]
        for row in accelerator['rows']:
            val=z(row['value']);reference=tuple(map(D,row['reference']))
            err=norm(sub(val,reference));rel=err/norm(reference)
            close(err,D(row['absolute_error']),'reference absolute error','1e-35')
            close(rel,D(row['relative_error']),'reference relative error','1e-30')
            if row['kind']=='influential':req(rel<D('1e-9'),'influential kernel accuracy')
            if row['kind']=='max_condition':bad.append(rel)
        req(len(bad)==1 and bad[0]>D('.01'),'retain ill-conditioned relative failure')
        for r in summary['rows']:
            raw=read(r['label']+'.json')
            for i,(value,saved) in enumerate(zip(r['values'],raw['values'])):
                err=norm(sub(z(value),z(saved)))
                req(err<D('2e-17'),'independent E/W replay')
                close(err,f(r['replay_absolute_difference'][i]),'replay residual','1e-20')
    return errors
