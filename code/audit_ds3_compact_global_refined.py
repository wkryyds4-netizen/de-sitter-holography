"""Fix structural SymPy equality in the original audit, retaining its failure."""
import hashlib,json
from pathlib import Path
import sympy as sp
from audit_ds3_compact_global import audit

if __name__=='__main__':
    root=Path(__file__).resolve().parents[1];base=root/'results/stage4/ds3-compact-global'
    dest=base/'audit-refined.json'
    if dest.exists():raise RuntimeError('Preserve audit')
    old=json.loads((base/'audit.json').read_text())
    result=audit(root)
    assert old['errors']==result['errors']==['central coefficient induction']
    n=sp.symbols('n',positive=True)
    lhs=(n+1)**3-(n+sp.Rational(1,2))**2*(n+2)
    rhs=(3*n+2)/4
    assert sp.expand(lhs-rhs)==0
    result['errors']=[];result['status']='passed';result['checks']+=2
    result['correction']='Original checker used structural == on factored versus distributed expressions. Exact expand(lhs-rhs)==0 verifies the unchanged induction identity; no data or tolerance changed.'
    result['hashes']={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in ['code/audit_ds3_compact_global.py','code/audit_ds3_compact_global_refined.py','results/stage4/ds3-compact-global/audit.json']}
    dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print(result['status'],result['checks'])
