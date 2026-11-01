"""Replay evidence integrity; deliberately does not mark pentagon as passed."""
from pathlib import Path
import json,hashlib,math
def audit(root):
    base=root/'results/stage4/unitary-flip';errors=[];count=0
    def ck(b,msg):
        nonlocal count
        count+=1
        if not b:errors.append(msg)
    study=json.loads((base/'study.json').read_text());diag=json.loads((base/'diagnostics.json').read_text())
    rows=[json.loads(s) for s in (base/'runs.jsonl').read_text().splitlines()]
    ck(rows==study['rows'],'raw runs disagree');ck(len(rows)==21,'run count')
    ck('unresolved' in study['status'],'missing scientific failure status')
    ck(len(diag['checks'])==36 and all(r['passed'] for r in diag['checks']),'scalar/algebra checks')
    for record in [study,diag]:
        for p,h in record['hashes'].items():ck(hashlib.sha256((root/p).read_bytes()).hexdigest()==h,'hash '+p)
    for r in rows:
        z=complex(*r['overlap']);n0,n1,n2,n3=r['norms']
        ck(all(abs(v-1)<2e-12 for v in [n0,n1,n2,n3]),'finite norm')
        ck(abs(r['D']-(n1*n1+n2*n2-2*z.real))<2e-12,'distance arithmetic')
        trace=n1*n1+n2*n2;disc=math.hypot(n1*n1-n2*n2,2*abs(z))
        expected=[(trace-disc)/2,(trace+disc)/2]
        ck(max(abs(a-b) for a,b in zip(r['gram_eigenvalues'],expected))<2e-12,'Gram arithmetic')
        ck(expected[0]>=-2e-12,'finite Gram positivity')
        ck(abs(complex(*r['projected_residual']))<=r['pentagon_residual']+2e-12,'projection bound')
        ck(all(t['phase_tail_bound']<1e-14 for t in r['product_tails']),'product tail')
        ck('unresolved' in r['residual_status'],'unresolved row status')
    # These identify retained observations, not acceptance criteria for a theorem.
    for method,n in [('periodic_fft',1024),('hermite_dvr',512)]:
        records=[r for r in rows if r['method']==method and r['n']==n]
        ck(len(records)==3,'terminal grid count')
        ck(next(r for r in records if r['flat'])['pentagon_residual']<1e-10,'flat control')
        ck(all(r['pentagon_residual']>.02 for r in records if not r['flat']),'non-flat residual unexpectedly lost')
    return dict(status='passed' if not errors else 'failed',checks=count,errors=errors,
        scientific_status='non-flat pentagon reproduction unresolved',
        scope='Hash and arithmetic integrity only; no continuum error bound or theorem refutation.')
if __name__=='__main__':
    root=Path(__file__).resolve().parents[1];p=root/'results/stage4/unitary-flip/audit.json'
    if p.exists():raise RuntimeError('Do not overwrite')
    r=audit(root);r['hashes']={'code/audit_unitary_flip.py':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    p.write_text(json.dumps(r,indent=2));print(json.dumps(r))
