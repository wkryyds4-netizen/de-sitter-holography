"""Verify the energy-weighted charged four-point reconstruction from saved samples."""
from pathlib import Path
import hashlib
import json
import numpy as np
from complex_syk_dictionary import fermions

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-dictionary'


def verify(root=ROOT):
    dest=root/'results/susskind-dictionary'
    record=json.loads((dest/'calibration.json').read_text(encoding='utf-8'))
    rows=[]
    for sample in record['samples']:
        n,p=sample['N'],sample['p']
        cs=fermions(n)
        with np.load(dest/sample['file'],allow_pickle=False) as data:
            h,a=data['H'],data['A2']
        sectors=[]
        for m in range(n+1):
            inds=np.array([s for s in range(1<<n) if s.bit_count()==m])
            e,u=np.linalg.eigh(h[np.ix_(inds,inds)])
            sectors.append((inds,e,u))
        dynamic=0.
        error=0.
        for m in range(1,n+1):
            lo,el,ul=sectors[m-1]
            hi,eh,uh=sectors[m]
            t=np.zeros((len(hi),len(hi)),dtype=complex)
            for c in cs:
                c_e=ul.conj().T@c[np.ix_(lo,hi)]@uh
                t+=c_e.conj().T@((el**2)[:,None]*c_e)
            direct=uh.conj().T@a[np.ix_(hi,hi)]@uh
            reconstructed=t-(m-2*p)*np.diag(eh**2)
            error=max(error,float(np.linalg.norm(direct-reconstructed)/max(1.,np.linalg.norm(direct))))
            gaps=eh[None,:]-eh[:,None]
            dynamic+=float(np.sum(np.abs(t[np.abs(gaps)>sample['zero_cutoff']])**2)/(1<<n))
        spectral_error=abs(dynamic-sample['dynamic_weight'])/max(1.,sample['G'])
        rows.append(dict(file=sample['file'],matrix_error=error,
                         fourpoint_dynamic_weight=dynamic,saved_dynamic_weight=sample['dynamic_weight'],
                         spectral_error=spectral_error,passed=error<2e-11 and spectral_error<2e-11))
    return dict(passed=all(r['passed'] for r in rows),rows=rows,
                scope='Coherent intermediate-state sum in same-sample neighboring charge sectors; finite-system check of equations (12)-(15).')


if __name__=='__main__':
    result=verify()
    result['code_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    with (OUT/'charge-bridge.json').open('x',encoding='utf-8') as f:
        json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({'passed':result['passed'],'samples':len(result['rows']),
                      'max_matrix_error':max(r['matrix_error'] for r in result['rows']),
                      'max_spectral_error':max(r['spectral_error'] for r in result['rows'])}))
    if not result['passed']:
        raise SystemExit(1)
