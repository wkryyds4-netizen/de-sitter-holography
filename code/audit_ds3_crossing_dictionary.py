"""Decimal replay of saved complex normalization, errors and tail bounds.

Does not rerun special functions or certify the Gauss remainder.
"""
import ast
import hashlib
import json
from decimal import Decimal as D, localcontext

PI=D('3.1415926535897932384626433832795028841971693993751058209749445923078164')


def audit(root):
    base=root/'results/stage4/ds3-crossing-dictionary'
    read=lambda name:json.loads((base/name).read_text(encoding='utf-8'))
    main,geo,orig,refined=map(read,['diagnostics.json','geometry.json','euclidean.json','euclidean-refined.json'])
    fix=read('source-correction.json');errors=[]
    def req(ok,msg):
        if not ok:errors.append('crossing dictionary: '+msg)
    def z(x):return (D(x['re']),D(x['im'])) if isinstance(x,dict) else (D(x),D(0))
    def add(a,b):return a[0]+b[0],a[1]+b[1]
    def sub(a,b):return a[0]-b[0],a[1]-b[1]
    def mul(a,b):return a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]
    def norm(a):return (a[0]**2+a[1]**2).sqrt()
    def scale(a,b):return a[0]*b,a[1]*b
    def div(a,b):return scale(mul(a,(b[0],-b[1])),1/(b[0]**2+b[1]**2))
    def close(a,b,msg,tol='1e-28'):
        req(norm(sub(a,b))<=D(tol)*max(norm(b),D(1)),msg)
    def sh(x):return (x.exp()-(-x).exp())/2
    def ch(x):return (x.exp()+(-x).exp())/2
    def K(a,b,p):
        denom=(ch(PI*(a+b))+ch(PI*p))*(ch(PI*(a-b))+ch(PI*p))/4
        return PI*sh(PI*a)*sh(PI*b)*sh(PI*p)/(8*a*b*p*denom)
    def stripped(text):
        tree=ast.parse(text)
        for node in ast.walk(tree):
            if isinstance(node,(ast.Module,ast.FunctionDef,ast.ClassDef,ast.AsyncFunctionDef)) and node.body:
                first=node.body[0]
                if isinstance(first,ast.Expr) and isinstance(first.value,ast.Constant) and isinstance(first.value.value,str):node.body=node.body[1:]
        return ast.dump(tree,include_attributes=False)
    old=root/fix['executed_source'];new=root/fix['current_source']
    req(hashlib.sha256(old.read_bytes()).hexdigest()==fix['executed_sha256'],'executed source hash')
    req(hashlib.sha256(new.read_bytes()).hexdigest()==fix['current_sha256'],'corrected source hash')
    req(stripped(old.read_text(encoding='utf-8'))==stripped(new.read_text(encoding='utf-8')),'executable AST unchanged')
    for data in [main,refined,geo]:
        for path,sha in data['hashes'].items():
            actual=root/path
            if path==fix['current_source'] and sha==fix['executed_sha256']:actual=old
            req(hashlib.sha256(actual.read_bytes()).hexdigest()==sha,'runtime hash '+path)
    req(main['passed']==main['total']==len(main['checks'])==33,'normalization check count')
    req(geo['passed']==geo['total']==len(geo['checks'])==10,'geometry check count')
    req(refined['status']=='completed' and len(refined['rows'])==6 and len(orig['rows'])==4,'complete quadrature grids')
    with localcontext() as ctx:
        ctx.prec=70
        for c in main['checks']:
            err=norm(sub(z(c['lhs']),z(c['rhs'])))/max(norm(z(c['rhs'])),D('1e-40'))
            req(c['passed'] is True and err<D(c['tolerance']),'saved check '+c['name'])
            close(z(c['error']),z(err),'recorded normalized residual','1e-40')
        for r in main['rows']:
            a,q,c,p=map(D,r['external']);b,d=map(D,[r['b'],r['d']])
            ws=2*b*b/PI*K(a,b,q)*K(b,c,p)
            wt=2*d*d/PI*K(c,d,q)*K(d,a,p)
            close(z(ws),z(r['ws']),'s rigging weight')
            close(z(wt),z(r['wt']),'t rigging weight')
            ns,nt=z(r['ns']),z(r['nt'])
            close(z(norm(ns)**2/ws),z(32*PI**7/b**2),'s raw norm')
            close(z(norm(nt)**2/wt),z(32*PI**7/d**2),'t raw norm')
            phase=div(scale(ns,1/norm(ns)),scale(nt,1/norm(nt)))
            u=scale(mul(phase,z(r['raw'])),b*d/PI**4)
            close(u,z(r['u']),'physical complex U')
            reverse=z(r['reverse']);close(u,(reverse[0],-reverse[1]),'Hermitian exchange')
        for r in orig['rows']+refined['rows']:
            analytic=z(r.get('analytic',orig['analytic']))
            value=z(r['value']);total=(D(0),D(0))
            req(len(r['charts'])==6,'six raw charts')
            for chart in r['charts']:total=add(total,z(chart))
            close(total,value,'sum of complex charts')
            err=norm(sub(value,analytic))/norm(analytic)
            close(z(err),z(r['relative_error']),'complex Euclidean residual')
            close(z(D(r['absolute'])/norm(value)),z(r['cancellation']),'cancellation ratio')
            req(D(r['absolute'])>=norm(value),'absolute sum lower bound')
        req(D(refined['rows'][1]['relative_error'])<D('1e-9'),'first refined complex check')
        req(D(refined['rows'][5]['relative_error'])<D('1e-7'),'second refined complex check')
        for r in geo['areas']:
            close(z(PI/4),z(r['expected']),'independent area target')
            err=abs(D(r['value'])/(PI/4)-1)
            req(err<D('1e-14'),'chart area coverage')
        for r in geo['tails']:
            ep=D(r['epsilon']);delta=ep/(1-ep);v=12-(delta.ln())
            bound=128*PI**3/3*(2*delta*(v*v+2*v+2)+delta*delta/2*(v*v+v+D('.5')))
            close(z(bound),z(r['bound']),'tail expression')
            close(z(bound),z(r['direct']),'tail primitive quadrature')
    return errors


if __name__=='__main__':
    from pathlib import Path
    errors=audit(Path(__file__).resolve().parents[1])
    print(json.dumps(dict(status='failed' if errors else 'passed',errors=errors),indent=2))
    if errors:raise SystemExit(1)
