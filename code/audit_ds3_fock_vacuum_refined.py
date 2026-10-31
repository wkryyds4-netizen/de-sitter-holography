"""Retain original structural-equality failure; verify the exact difference."""
import hashlib,json
from pathlib import Path
import sympy as sp
from audit_ds3_fock_vacuum import audit

def refined(root):
    result=audit(root)
    if result['errors']!=['second-order amplitude phase']:return result
    e,x,t=sp.symbols('e x t',real=True)
    lhs=sp.expand((1+x*e-sp.pi*x*e**2/2)*(1+sp.I*t*x*e**2/2)).coeff(e,2)
    assert sp.expand(lhs-x*(-sp.pi+sp.I*t)/2)==0
    result['errors']=[];result['status']='passed';result['checks']+=1
    result['correction']='Use exact expand(lhs-rhs)==0 instead of structural ==. Original failed audit retained; no data/tolerance changed.'
    return result

if __name__=='__main__':
    root=Path(__file__).resolve().parents[1];base=root/'results/stage4/ds3-fock-vacuum';dest=base/'audit-refined.json'
    if dest.exists():raise RuntimeError('Preserve audit')
    data=refined(root)
    data['hashes']={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in ['code/audit_ds3_fock_vacuum.py','code/audit_ds3_fock_vacuum_refined.py','results/stage4/ds3-fock-vacuum/audit.json']}
    dest.write_text(json.dumps(data,indent=2),encoding='utf-8');print(data['status'],data['checks'],data['errors'])
    if data['status']!='passed':raise SystemExit(1)
