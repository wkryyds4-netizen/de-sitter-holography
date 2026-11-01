"""Replay saved Galerkin coefficients and provenance, not a continuum theorem."""
from pathlib import Path
import json,hashlib,math
import numpy as np

def audit(root):
    base=root/'results/stage4/flip-domain';errors=[];count=0
    def ck(ok,msg):
        nonlocal count
        count+=1
        if not ok:errors.append(msg)
    study=json.loads((base/'study.json').read_text())
    raw=[json.loads(line) for line in (base/'runs.jsonl').read_text().splitlines()]
    records=[study]+[json.loads((base/(n+'.json')).read_text()) for n in ['diagnostics','refinement','regularity']]
    ck(raw==study['rows'],'raw/study equality');ck(len(raw)==11,'11 runs')
    ck('pending' in study['status'],'uncertified status')
    for record in records:
        for p,h in record['hashes'].items():
            ck(hashlib.sha256((root/p).read_bytes()).hexdigest()==h,'hash '+p)
    diag,ref,reg=records[1:]
    ck(len(diag['checks'])==22,'initial check count')
    ck(sum(not c['passed'] for c in diag['checks'])==2 and diag['status']=='failed','preserved initial failures')
    ck(len(ref['checks'])==5 and all(c['passed'] for c in ref['checks']),'refinement')
    ck(len(reg['checks'])==31 and all(c['passed'] for c in reg['checks']),'regularity controls')
    for i,row in enumerate(raw):
        p=root/row['coefficients'];ck(hashlib.sha256(p.read_bytes()).hexdigest()==row['coefficients_sha256'],f'npz hash {i}')
        ck(row['q']>row['n'],f'oversampling {i}')
        with np.load(p,allow_pickle=False) as saved:
            branches={key:saved[key] for key in ['left','right']}
        all_leaks=[]
        for key,length in [('left',3),('right',6)]:
            vec=branches[key];n=row['n'];ck(vec.shape==(length,n,n),f'shape {i} {key}')
            vacuum=np.zeros((n,n),complex);vacuum[0,0]=1
            ck(np.array_equal(vec[0],vacuum),f'input {i} {key}')
            norms=[math.sqrt(math.fsum(float(z.real)**2+float(z.imag)**2 for z in v.flat)) for v in vec]
            for j,(got,want) in enumerate(zip(norms,row[key+'_norms'])):
                ck(abs(got-want)<3e-13,f'norm {i} {key} {j}')
            leaks=[math.sqrt(max(0,a*a-b*b)) for a,b in zip(norms[:-1],norms[1:])]
            ck(all(b<=a+3e-13 for a,b in zip(norms[:-1],norms[1:])),f'contraction {i} {key}')
            for j,(got,want) in enumerate(zip(leaks,row[key+'_leaks'])):
                # Near-zero loss: propagate 6e-13 squared-norm allowance through sqrt.
                ck(abs(got*got-want*want)<6e-13,f'leak squared {i} {key} {j}')
            all_leaks+=leaks
        delta=branches['left'][-1]-branches['right'][-1]
        residual=math.sqrt(math.fsum(float(z.real)**2+float(z.imag)**2 for z in delta.flat))
        ck(abs(residual-row['residual'])<3e-13,f'residual {i}')
        ck(abs(delta[0,0]-complex(*row['vacuum_projection']))<3e-13,f'vacuum projection {i}')
        ck(abs(residual-math.fsum(all_leaks)-row['quadrature_dependent_margin'])<7*math.sqrt(6e-13),f'margin {i}')
        ck('not a continuum lower bound' in row['scope'],f'scope {i}')
    return dict(status='passed' if not errors else 'failed',checks=count,errors=errors,
                scientific_status='bounded pentagon unresolved; quadrature not certified',
                scope='Saved-vector arithmetic and hash replay. Does not recalculate matrix elements or certify continuum claims.')

if __name__=='__main__':
    root=Path(__file__).resolve().parents[1];out=root/'results/stage4/flip-domain/audit.json'
    if out.exists():raise RuntimeError('Do not overwrite')
    data=audit(root);data['hashes']={'code/audit_flip_domain.py':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    out.write_text(json.dumps(data,indent=2));print(json.dumps(data))
