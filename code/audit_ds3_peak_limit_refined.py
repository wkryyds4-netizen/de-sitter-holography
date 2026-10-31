"""Repair audit precision using an independently summed Machin formula for pi.

The initial audit and its failed output remain unchanged. Scientific data
and tolerances are unchanged; only the insufficient 70-digit pi is replaced.
"""
import hashlib,json
from pathlib import Path
from decimal import Decimal as D,localcontext
import audit_ds3_peak_limit as original


def machin():
    with localcontext() as ctx:
        ctx.prec=180
        def atan_inverse(n):
            x=1/D(n);power=x;total=D(0);k=0
            while True:
                term=power/(2*k+1)
                total+=term if k%2==0 else -term
                power*=x*x;k+=1
                next_term=power/(2*k+1)
                if next_term<D('1e-165'):return total,next_term,k
        a,ea,na=atan_inverse(5);b,eb,nb=atan_inverse(239)
        return 16*a-4*b,16*ea+4*eb,dict(atan5_terms=na,atan239_terms=nb)


if __name__=='__main__':
    root=Path(__file__).resolve().parents[1];dest=root/'results/stage4/ds3-peak-limit/audit-refined.json'
    if dest.exists():raise RuntimeError('Preserve prior refined audit')
    old_pi=original.PI;pi,bound,terms=machin();original.PI=pi
    assert abs(pi-old_pi)<D('1e-70') and bound<D('1e-160')
    result=original.audit(root)
    result['pi_refinement']=dict(pi=str(pi),series_truncation_bound=str(bound),terms=terms,
          reason='Initial imported pi had only about 70 decimal places; the unchanged audit demanded 85-90 decimal places.',
          scope='Alternating-series truncation bound only; Decimal rounding separately limited by 180-digit evaluation.')
    paths=['code/audit_ds3_peak_limit.py','code/audit_ds3_peak_limit_refined.py','results/stage4/ds3-peak-limit/audit.json']
    result['hashes']={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in paths}
    dest.write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2));raise SystemExit(result['status']!='passed')
