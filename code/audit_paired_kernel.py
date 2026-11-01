"""Read-only replay of the paired-kernel dictionary and saved heat arithmetic."""
from pathlib import Path
import json,hashlib
import mpmath as m

def audit(root):
    base=root/'results/stage4/paired-kernel';count=0;errors=[]
    def ck(ok,msg):
        nonlocal count
        count+=1
        if not ok:errors.append(msg)
    study=json.loads((base/'heat-study.json').read_text());diag=json.loads((base/'diagnostics.json').read_text())
    raw=[json.loads(x) for x in (base/'heat-runs.jsonl').read_text().splitlines()]
    ck(raw==study['rows'],'raw/study equality');ck(len(raw)==240,'240 nested-prefix rows')
    ck(len({(r['seed'],r['h'],r['delta'],r['power']) for r in raw})==240,'unique sampling settings')
    ck(len(diag['checks'])==30 and all(c['passed'] for c in diag['checks']),'30 dictionary/control checks')
    for d in [study,diag]:
        for p,h in d['hashes'].items():ck(hashlib.sha256((root/p).read_bytes()).hexdigest()==h,'hash '+p)
    with m.workdps(70):
        for i,r in enumerate(raw):
            t=m.mpf(str(r['t']));delta=m.mpf(str(r['delta']));mean=m.mpf(str(r['mean_phase_loss']));v2=m.mpf(str(r['mean_square']))
            ck(abs(t-2*delta)<m.mpf('1e-18') and r['n']==2**r['power'],f'parameter relation {i}')
            ck(0<=mean<=2 and mean**2-m.mpf('1e-14')<=v2<=2*mean+m.mpf('1e-14'),f'bounded moments {i}')
            ck(0<=r['count_above_point1']<=r['n'] and 0<=r['maximum']<=2.000000000001,f'loss range {i}')
            heat=(m.pi*t+mean)/(1+m.pi*t)
            error=2*m.sqrt(2*(m.pi*delta)**2/((1+m.pi*delta)*(1+2*m.pi*delta)))+m.sqrt(heat)
            ck(abs(heat-m.mpf(str(r['estimated_one_minus_C'])))<m.mpf('1e-14'),f'heat arithmetic {i}')
            ck(abs(error-m.mpf(str(r['estimated_weak_error_functional'])))<m.mpf('1e-14'),f'error-functional arithmetic {i}')
            if r['h'] is None:
                exact=2*m.pi*t/(1+2*m.pi*t)
                ck(abs(exact-m.mpf(str(r['exact_one_minus_C'])))<m.mpf('1e-14'),f'chirp expectation {i}')
        dec=lambda z:m.mpc(*z)
        for i,r in enumerate(diag['scalar_rows']):
            ck(abs(abs(dec(r['direct'])-dec(r['factored']))-m.mpf(r['factor_error']))<m.mpf('1e-48'),f'factor replay {i}')
            ck(abs(abs(dec(r['inversion'])-dec(r['target']))-m.mpf(r['inversion_error']))<m.mpf('1e-48'),f'inversion replay {i}')
        for i,r in enumerate(diag['noninteger_rows']):
            ck(abs(abs(dec(r['value'])-dec(r['branch_target']))-m.mpf(r['difference']))<m.mpf('1e-38'),f'noninteger replay {i}')
    return dict(status='passed' if not errors else 'failed',checks=count,errors=errors,
        scientific_status='bounded-kernel formulation and empirical heat-tail study; pentagon not evaluated',
        scope='Saved arithmetic and hash replay, not a rerun of Sobol integration or certification of its error.')

if __name__=='__main__':
    root=Path(__file__).resolve().parents[1];p=root/'results/stage4/paired-kernel/audit.json'
    if p.exists():raise RuntimeError('Do not overwrite')
    d=audit(root);d['hashes']={'code/audit_paired_kernel.py':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    p.write_text(json.dumps(d,indent=2));print(json.dumps(d))
