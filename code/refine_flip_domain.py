"""Keep initial failures; check exact difference and a smaller pole shell."""
from pathlib import Path
import json,hashlib
import sympy as s
import mpmath as m
root=Path(__file__).resolve().parents[1];base=root/'results/stage4/flip-domain'
out=base/'refinement.json'
if out.exists():raise RuntimeError('Do not overwrite refinement')
original=json.loads((base/'diagnostics.json').read_text());checks=[]
def ck(name,ok):checks.append(dict(name=name,passed=bool(ok)))
ck('preserved two initial failures',sum(not r['passed'] for r in original['checks'])==2)
a,b=s.symbols('a b',positive=True);x=s.pi*(a*a+a+b*b)/b;y=-s.pi*(a+1)/b
diff=s.simplify(x*y+s.pi**2+s.pi**2*a*((a+1)**2+b*b)/b**2)
ck('wedge identity by exact difference',diff==0)
m.mp.dps=60;h=m.mpf('.5');q=m.exp(-m.pi*h);qt=m.exp(-m.pi/h)
R=-m.qp(qt*qt,qt*qt)/m.qp(q*q,q*q)
C=2/(m.pi*h)*abs(R)**2*m.exp(-m.pi*(h*h+h**-2));rows=[]
for nodes in [64,128]:
    r=m.mpf('.0125');total=0
    for j in range(nodes):
        t=2*m.pi*j/nodes;z=m.pi*h-1j*m.pi+r*m.exp(1j*t)
        ph=m.qp(-qt*m.exp(z/(1j*h)),qt*qt)/m.qp(-q*m.exp(z),q*q)
        xx=m.pi*h+r*m.cos(t);yy=-m.pi/h+r*m.sin(t)/h
        total+=abs(ph)**2*m.exp(-(xx*xx+yy*yy)/m.pi)/m.pi**2
    shell=r*r/h*2*m.pi/nodes*total
    rows.append(dict(nodes=nodes,shell=m.nstr(shell,45),relative_to_limit=m.nstr(shell/C-1,40)))
    ck('same one-percent criterion on smaller shell '+str(nodes),abs(shell/C-1)<m.mpf('.01'))
ck('angular refinement agreement',abs(m.mpf(rows[0]['shell'])/m.mpf(rows[1]['shell'])-1)<m.mpf('1e-35'))
result=dict(status='passed' if all(r['passed'] for r in checks) else 'failed',checks=checks,
    initial_status=original['status'],reason='SymPy structural equality replaced by zero difference; h=.5 radius halved, original 1% criterion unchanged.',
    h='.5',radius='.0125',precision=60,rows=rows,
    scope='Repairs two finite checks; does not certify Galerkin quadrature or refute the pentagon.',
    hashes={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in ['code/refine_flip_domain.py','results/stage4/flip-domain/diagnostics.json']})
out.write_text(json.dumps(result,indent=2));print(result['status'],len(checks))
