"""Symbolic degree-two Gram audit from the two-pair limiting rule.

Reads sealed finite CAR data, never imports/reruns its generator.
The symbolic audit is algebraic consistency, not an all-word theorem.
"""
from pathlib import Path
import hashlib, json, platform
from datetime import datetime, timezone
from itertools import product
import sympy as S

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-slater-joint-fourth/gram-checks.json'
SOURCE=ROOT/'results/susskind-slater-joint-fourth/exact-checks.json'
t=S.Symbol('gamma', positive=True)
q=t**4; s=t**2
PAIRINGS=[((0,1),(2,3)),((0,2),(1,3)),((0,3),(1,2))]
checks=[]
def check(name,lhs,rhs):
    residual=S.factor(lhs-rhs)
    assert residual==0,(name,residual)
    checks.append(name)
def kind(w,e):
    a,b=(w[i] for i in e)
    if a==b=='H':return 'G'
    if a==b and a in '+-':return 'X'
    return 'B' if a=='+' or b=='-' else 'S'
def moment(w):
    if not w:return S.Integer(1)
    if len(w)%2:return S.Integer(0)
    if len(w)==2:return {'HH':1,'H+':t,'-H':t,'-+':1}.get(w,S.Integer(0))
    assert len(w)==4
    total=0
    for pi,pairs in enumerate(PAIRINGS):
        kinds=[kind(w,e) for e in pairs]
        if any(k in ('B','X') for k in kinds):continue
        v=t**sum(w[i]!=w[j] and 'H' in (w[i],w[j]) for i,j in pairs)
        if pi==1:v*=q
        if pi==2 and kinds[1]=='S':v*=q
        total+=v
    return S.expand(total)
def add(*polys):
    out={}
    for poly in polys:
        for w,c in poly.items():out[w]=S.expand(out.get(w,0)+c)
    return {w:c for w,c in out.items() if c!=0}
def scale(c,p):return {w:S.expand(c*v) for w,v in p.items()}
def mul(a,b):
    return add(*({w+v:c*d} for w,c in a.items() for v,d in b.items()))
def adj(w):return ''.join({'+':'-','-':'+','H':'H'}[x] for x in reversed(w))
def inner(a,b,fn=moment):return S.expand(sum(c*d*fn(adj(w)+v) for w,c in a.items() for v,d in b.items()))
I={'':1}; H={'H':1}; P={'+':1}; A={'-':1}
U={'++':1}
V=add({'+H':1},scale(-t,U))
Dv={'H+':1,'+H':-1,'':-t}
Rv=add({'HH':1,'':-1},scale(-s,U),scale(-2*t,V),scale(-t,Dv))
vectors=[I,P,add(H,scale(-t,P)),U,V,Dv,Rv]
norms=[1,1,1-s,2*q,q*(1-s),1-q,(1-s)**2*(1+s+2*q)]
for i,a in enumerate(vectors):
    for j,b in enumerate(vectors):
        check(f'Gram {i},{j}',inner(a,b),norms[i] if i==j else 0)
# Residual vectors have degree <=2, so their norms are determined here.
def null(name,a):check(name,inner(a,a),0)
null('he1 action',add(mul(H,P),scale(-t,I),scale(-t,U),scale(-1,V),scale(-1,Dv)))
null('h residual action',add(mul(H,add(H,scale(-t,P))),scale(-(1-s),I),scale(-t,V),scale(-1,Rv)))
null('Mh vacuum',add(mul(A,H),scale(-t,I)))
g=add(H,scale(-t,add(P,A)))
g2=mul(g,g)
check('g variance',inner(g,g),1-s)
check('g fourth',inner(g2,g2),(1-s)**2*(2+s+2*q))
null('centered g2 is Rv',add(g2,scale(-(1-s),I),scale(-1,Rv)))
check('h3 Mdag',moment('HHH+'),t*(2+q))
check('h2 Mdag2',moment('HH++'),2*s*q)
check('h Mdag h Mdag',moment('H+H+'),s)
finite=json.loads(SOURCE.read_text(encoding='utf-8'))
finite_checks=[]
for fixture in finite['fixtures']:
    L,D=fixture['L'],fixture['D']
    raw={row['word']:row['total'] for row in fixture['words']}
    raw.update(fixture['second_moments']);raw['']=1
    def fm(w):
        if len(w)%2:return S.Integer(0)
        return S.Rational(raw[w])/S.sqrt(S.Integer(D)**w.count('H')*S.Integer(L)**(len(w)-w.count('H')))
    v={'-H':1,'':-S.sqrt(S.Rational(L,D))}
    check(f'finite Mh residual {fixture["N"]},{fixture["p"]}',inner(v,v,fm),S.Rational(1,D))
    finite_checks.append(dict(N=fixture['N'],p=fixture['p'],squared_residual=str(S.Rational(1,D))))
if OUT.exists():raise FileExistsError(OUT)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
result=dict(status='passed',created_utc=datetime.now(timezone.utc).isoformat(),
    parameters=dict(q='gamma^4',s='gamma^2',domain='0 < gamma < 1',max_word_degree=4),
    checks=checks,count=len(checks),vectors=vectors,norms=list(map(str,norms)),finite_checks=finite_checks,
    environment=dict(python=platform.python_version(),sympy=S.__version__),
    input_sha256={p.relative_to(ROOT).as_posix():sha(p) for p in (Path(__file__),SOURCE)},
    scope='Exact symbolic consequences of stated two-pair limit; four finite residuals read from independent sealed CAR data. No higher-word or dynamical closure claim.')
with OUT.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2,default=str);f.write('\n')
print(json.dumps(dict(status='passed',checks=len(checks),sha256=sha(OUT))))
