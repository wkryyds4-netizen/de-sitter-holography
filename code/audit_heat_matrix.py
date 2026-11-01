"""Exact certificate replay and saved matrix arithmetic; not quadrature proof."""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib
import mpmath as m
from rational_heat_certificate import certify

def audit(root):
    base=root/'results/stage4/heat-matrix';errors=[];count=0
    def ck(ok,msg):
        nonlocal count
        count+=1
        if not ok:errors.append(msg)
    cert=json.loads((base/'rational-certificate.json').read_text())
    iv=json.loads((base/'cutoff-bound.json').read_text())
    study=json.loads((base/'matrix-study.json').read_text());direct=json.loads((base/'direct-left.json').read_text())
    raw=[json.loads(x) for x in (base/'matrix-runs.jsonl').read_text().splitlines()]
    ck(raw==study['rows'],'raw/study equality');ck(len(raw)==37,'37 matrix records');ck(len(direct['rows'])==12,'12 direct-left records')
    for data in [cert,iv,study,direct]:
        for p,h in data['hashes'].items():ck(hashlib.sha256((root/p).read_bytes()).hexdigest()==h,'hash '+p)
    def dyadic(t):
        sign,mant,exp,bits=t
        ck(mant>=0 and (mant.bit_length()==bits or mant==0),'dyadic mantissa')
        return (-1 if sign else 1)*F(mant)*F(2)**exp
    for row in cert['rows']:
        replay=certify(row['epsilon'])
        ck(replay==row,'exact rational certificate replay '+row['epsilon'])
        enclosure=row['values']['total'];lo=F(*enclosure['lower']);hi=F(*enclosure['upper'])
        ck(lo<=hi<F(*row['strict_rational_upper']),'strict rational upper '+row['epsilon'])
        for a in [x for x in iv['rows'] if x['epsilon']==row['epsilon']]:
            pair=a['values']['total_matrix_element_cutoff_bound']['dyadic_endpoints']
            l,h=map(dyadic,pair)
            ck(lo<=l<=h<=hi,'independent interval enclosure '+str(a['dps'])+row['epsilon'])
    with m.workdps(70):
        dec=lambda z:m.mpc(str(z[0]),str(z[1]))
        for i,r in enumerate(raw):
            defect=dec(r['lhs'])-dec(r['rhs'])
            ck(abs(defect-dec(r['defect']))<m.mpf('3e-16'),f'complex subtraction {i}')
            ck(abs(abs(defect)-m.mpf(str(r['defect_modulus'])))<m.mpf('3e-16'),f'modulus {i}')
            ck(abs(defect)<=m.mpf(str(r['word_vector_residual']))+m.mpf('1e-12'),f'normalized projection {i}')
            ck(0<=r['left_norm']<=1.00000000001 and 0<=r['right_norm']<=1.00000000001,f'contractions {i}')
        for i,r in enumerate(direct['rows']):
            if r['unit']:
                expected=1/(1+m.pi*m.mpf(str(r['epsilon'])))
                ck(abs(dec(r['value'])-expected)<m.mpf('1e-12'),f'exact no-phase Gaussian {i}')
            if r['flat']:
                ref=next(x for x in raw if x['flat'] and x['method']=='periodic_fft' and x['epsilon']==r['epsilon'])
                ck(abs(dec(r['value'])-dec(ref['lhs']))<m.mpf('1e-12'),f'independent flat left control {i}')
    flat=[r for r in raw if r['flat']]
    ck(all(r['defect_modulus']<1e-12 for r in flat if r['epsilon']==0),'undamped flat control')
    ck(all(r['defect_modulus']>.03 for r in flat if r['epsilon']==.01),'damped flat identity not imposed')
    return dict(status='passed' if not errors else 'failed',checks=count,errors=errors,
        scientific_status='heat cutoff certified; matrix quadrature uncertified; no continuum counterexample',
        scope='Re-evaluates exact elementary rational certificate and saved complex arithmetic; does not certify FFT/Hermite/direct quadrature.')

if __name__=='__main__':
    root=Path(__file__).resolve().parents[1];out=root/'results/stage4/heat-matrix/audit.json'
    if out.exists():raise RuntimeError('Do not overwrite')
    d=audit(root);d['hashes']={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in
        ['code/audit_heat_matrix.py','code/rational_heat_certificate.py']}
    out.write_text(json.dumps(d,indent=2));print(json.dumps(d))
