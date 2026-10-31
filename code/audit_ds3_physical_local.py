"""Independent Decimal replay of physical weights on saved full-kernel nodes.

Does not import the implementation, mpmath, or sympy. Does not certify U itself.
"""
import hashlib,itertools,json
from decimal import Decimal as D,localcontext
from pathlib import Path


def audit(root):
    with localcontext() as ctx:
        ctx.prec=95
        errors=[];count=0
        def check(name,yes):
            nonlocal count
            count+=1
            if not yes:errors.append(name)
        def atan(x):
            s=D(0);term=x;n=0
            while abs(term)>D('1e-100'):
                s+=term/D(2*n+1);term*=-x*x;n+=1
            return s
        pi=16*atan(D(1)/5)-4*atan(D(1)/239)
        cosh=lambda x:(x.exp()+(-x).exp())/2
        h=lambda x:(1+x*x).sqrt()
        def g(x):return x*(pi/(h(x)+x)).exp()/(1+(-2*pi*x).exp())
        def bal(x,y,p):
            A=cosh(pi*(x+y));B=cosh(pi*(x-y));C=cosh(pi*p)
            return (A-B)/(8*pi*(A+C)*(B+C))
        def amp(M,u,v,p,q,r):
            a,b,c,d=M+u+v,M+u,M,M+v+r
            return g(a)/g(c)*(bal(a,b,q)*bal(b,c,p)*bal(a,d,p)*bal(d,c,q)).sqrt()
        def cis(x):
            co=D(1);si=x;ct=D(1);st=x;n=1
            while max(abs(ct),abs(st))>D('1e-100'):
                ct*=-x*x/D((2*n-1)*2*n);st*=-x*x/D(2*n*(2*n+1))
                co+=ct;si+=st;n+=1
            return co,si
        data=json.loads((root/'results/stage4/ds3-physical-local/diagnostics.json').read_text())
        check('driver successful',data['status']=='passed')
        check('expected record shapes',len(data['points'])==12 and len(data['derivatives'])==72 and len(data['integrals'])==12)
        for c in data['checks']:
            check('stored check '+c['name'],c['passed'] is True and D(str(c['error']))<=D(str(c['tolerance'])))
        for name,sha in data['hashes'].items():
            check('input hash '+name,hashlib.sha256((root/name).read_bytes()).hexdigest()==sha)
        mat=[[1,1,1,0],[1,1,0,0],[1,0,0,0],[1,0,1,1]]
        det=sum((-1)**sum(p[i]>p[j] for i in range(4) for j in range(i+1,4))*
                __import__('math').prod(mat[i][p[i]] for i in range(4)) for p in itertools.permutations(range(4)))
        check('signed determinant',det==data['jacobian_determinant']==-1)
        for i,row in enumerate(data['points']):
            M,u,v,p,q,r=[D(str(row[k])) for k in ['M','u','v','p','q','r']]
            actual=amp(M,u,v,p,q,r)
            check('Decimal amplitude '+str(i),abs(actual-D(row['amplitude']))<D('1e-70'))
            a,b,d=M+u+v,M+u,M+v+r
            check('energy Jacobian '+str(i),abs(h(a)*h(b)*h(d)/(a*b*d)-D(row['energy_measure_jacobian']))<D('1e-70'))
        nodes=[json.loads(x) for x in (root/'results/stage4/ds3-uniform-local/integral-nodes.jsonl').read_text().splitlines()]
        u,v,p,q=map(D,['.1','-.1','.4','.4']);replayed=[]
        for row in data['integrals']:
            M=D(row['M']);t=D(row['t']);sr=si=tail=D(0)
            selected=[n for n in nodes if n['M']==row['M'] and n['order']==row['order']]
            check('quadrature count '+str((M,t,row['order'])),len(selected)==row['order'])
            for n in selected:
                r,w,f=map(D,[n['r'],n['weight'],n['test']]);ur,ui=map(D,n['value'])
                A=amp(M,u,v,p,q,r)
                phase=t*(h(M+u+v)+h(M)-h(M+u)-h(M+v+r));cr,ci=cis(phase)
                scale=w*f*A;sr+=scale*(cr*ur-ci*ui);si+=scale*(cr*ui+ci*ur)
                tail+=abs(scale)*D(n['tail'])
            vr,vi=map(D,row['value'])
            check('Decimal integral '+str((M,t,row['order'])),max(abs(sr-vr),abs(si-vi))<D('1e-70'))
            check('series-tail propagation '+str((M,t,row['order'])),abs(tail-D(row['propagated_series_tail']))<D('1e-70')*max(tail,D('1e-300')))
            replayed.append(dict(M=int(M),order=row['order'],t=str(t),real=str(sr),imag=str(si)))
        refinements=[]
        for M in [64,256]:
            for t in ['0.0','1.0','2.0']:
                rows=[x for x in data['integrals'] if x['M']==M and x['t']==t]
                x,y=sorted(rows,key=lambda z:z['order'])
                dr=D(x['value'][0])-D(y['value'][0]);di=D(x['value'][1])-D(y['value'][1])
                diff=(dr*dr+di*di).sqrt()
                refinements.append(dict(M=M,t=t,absolute_difference=str(diff),relative_to_target=str(diff/D(y['target'][0]))))
        return dict(status='passed' if not errors else 'failed',checks=count,errors=errors,
                    decimal_precision=95,replayed_integrals=replayed,refinements=refinements,
                    scope='Independent elementary-weight replay on immutable old U nodes. Not a new kernel computation, interval arithmetic, or a global physical limit.')


if __name__=='__main__':
    root=Path(__file__).resolve().parents[1];dest=root/'results/stage4/ds3-physical-local/audit.json'
    if dest.exists():raise RuntimeError('Preserve old audit')
    result=audit(root);dest.write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps({k:result[k] for k in ['status','checks','errors']}))
    if result['errors']:raise SystemExit(1)
