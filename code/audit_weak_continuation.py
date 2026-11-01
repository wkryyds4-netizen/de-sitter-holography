"""Recompute saved complex arithmetic at 70 digits; no integrator is called."""
from pathlib import Path
import json,hashlib
import mpmath as m

def audit(root):
    base=root/'results/stage4/weak-continuation';errors=[];count=0
    def ck(ok,msg):
        nonlocal count
        count+=1
        if not ok:errors.append(msg)
    study=json.loads((base/'study.json').read_text());positive=json.loads((base/'positive-form.json').read_text())
    diag=json.loads((base/'diagnostics.json').read_text())
    raw=[json.loads(s) for s in (base/'runs.jsonl').read_text().splitlines()]
    ck(raw==study['rows'],'raw/study equality');ck(len(raw)==13,'13 records')
    ck(len(diag['checks'])==34 and all(c['passed'] for c in diag['checks']),'34 diagnostics')
    for record in [study,positive,diag]:
        for p,h in record['hashes'].items():ck(hashlib.sha256((root/p).read_bytes()).hexdigest()==h,'hash '+p)
    with m.workdps(70):
        dec=lambda z:m.mpc(*z);tol=m.mpf('1e-34');mass=m.pi*m.mpf(study['radius'])**2/(study['power']+1)
        ck(abs(m.mpf(positive['mass'])-mass)<tol,'radial mass')
        ck(abs(m.mpf(positive['norm'])-1)<tol,'positive norm')
        for i,row in enumerate(raw):
            if row['kind']=='integral':
                coarse,fine,direct=[dec(row[k]) for k in ['subtracted_coarse','subtracted_fine','direct']]
                ck(abs(abs(fine-coarse)-m.mpf(row['quadrature_difference']))<tol,f'quadrature difference {i}')
                ck(abs(abs(fine-direct)-m.mpf(row['direct_difference']))<tol,f'direct difference {i}')
            elif row['kind']=='wirtinger':
                v=list(map(dec,row['values']));delta=m.mpf(row['delta'])
                d=(v[0]-v[1]+1j*(v[2]-v[3]))/(4*delta)
                ck(abs(d-dec(row['derivative']))<tol,f'CR arithmetic {i}')
                ck(abs(dec(row['expected'])-m.pi)<tol,f'CR expectation {i}')
                ck(abs(abs(d/m.pi-1)-m.mpf(row['relative_error']))<m.mpf('1e-29'),f'CR error {i}')
                if 'direct_values' in row:
                    v=list(map(dec,row['direct_values']));dd=(v[0]-v[1]+1j*(v[2]-v[3]))/(4*delta)
                    ck(abs(dd-dec(row['direct_derivative']))<tol,f'direct CR {i}')
                    ck(abs(abs(dd-d)-m.mpf(row['direct_difference']))<tol,f'direct CR difference {i}')
            elif row['kind']=='morera':
                n=row['loop_nodes'];r=m.mpf(row['loop_radius']);v=list(map(dec,row['samples']))
                ck(len(v)==n,f'loop count {i}')
                result=m.fsum(v[k]*1j*r*m.exp(2j*m.pi*k/n)*2*m.pi/n for k in range(n))
                expected=2j*m.pi*mass
                ck(abs(result-dec(row['integral']))<tol,f'loop arithmetic {i}')
                ck(abs(expected-dec(row['expected']))<tol,f'loop expectation {i}')
                ck(abs(abs(result/expected-1)-m.mpf(row['relative_error']))<m.mpf('1e-29'),f'loop error {i}')
        expected=m.pi/(mass*(-2*m.pi+1j*m.pi))
        ck(len(positive['rows'])==3,'three positive-state differences')
        for i,row in enumerate(positive['rows']):
            v=list(map(dec,row['values']));delta=m.mpf(row['delta']);d=(v[0]-v[1]+1j*(v[2]-v[3]))/(4*delta)
            ck(abs(d-dec(row['derivative']))<m.mpf('1e-32'),f'positive derivative {i}')
            ck(abs(expected-dec(row['expected']))<m.mpf('1e-32'),f'positive expectation {i}')
            ck(abs(abs(d/expected-1)-m.mpf(row['relative_error']))<m.mpf('1e-29'),f'positive error {i}')
    return dict(status='passed' if not errors else 'failed',checks=count,errors=errors,
        scientific_status='literal local weak pairing is not holomorphic; bounded pentagon not tested',
        scope='70-digit replay of hashes and saved arithmetic; no certified quadrature or general operator-identity claim.')

if __name__=='__main__':
    root=Path(__file__).resolve().parents[1];out=root/'results/stage4/weak-continuation/audit.json'
    if out.exists():raise RuntimeError('Do not overwrite')
    data=audit(root);data['hashes']={'code/audit_weak_continuation.py':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    out.write_text(json.dumps(data,indent=2));print(json.dumps(data))
