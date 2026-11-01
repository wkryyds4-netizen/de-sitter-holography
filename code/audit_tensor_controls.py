"""Compare two enclosures by intersection, not enclosure-width ordering.

The initial control asked the narrower new enclosure to contain the old
one. Neither method requires that. Preserve that failed control verbatim.
"""
from pathlib import Path
import sys,json,hashlib
root=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(root/'tmp/arb-runtime'))
from flint import arb,fmpq,ctx
ctx.prec=256
def endpoints(item):
    return [arb(fmpq(int(a),int(b))) for a,b in [item['lower'],item['upper']]]
base=root/'results/stage4/residual-certificate'
old=json.loads((base/'tensor-controls.json').read_text())
pilot=json.loads((base/'rectangle-pilot.json').read_text())
checks=[dict(control=r['control'],passed=r['passed']) for r in old['checks'][:3]]
for part,row in zip(['real','imag'],old['checks'][3:]):
    lo,hi=endpoints(row['value']);error=endpoints(row['analytic_error'])[1]
    plo,phi=endpoints(pilot['finite_product_overlap_'+part])
    checks.append(dict(control='independent enclosures intersect: '+part,
                       passed=bool(lo-error<phi and hi+error>plo)))
output=base/'tensor-controls-audit.json'
payload=dict(passed=all(r['passed'] for r in checks),checks=checks,
             correction='Containment was an invalid acceptance criterion for independent enclosures of different widths; original failure preserved.',
             hashes={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in
             ['code/audit_tensor_controls.py','results/stage4/residual-certificate/tensor-controls.json',
              'results/stage4/residual-certificate/rectangle-pilot.json']})
with output.open('x') as f:json.dump(payload,f,indent=2)
print(json.dumps(payload))
assert payload['passed']
