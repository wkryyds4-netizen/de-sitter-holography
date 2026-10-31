"""Independent checks of saved pure-recoil results and limiting cases."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import mpmath as mp
from recoil_model import (recoil_hyper, recoil_series, recoil_angle,
    recoil_conformal, conformal_boundary, time_branch_q, q_from_time)

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'results'/'stage3'


def main():
    mp.mp.dps = 70
    study = json.loads((OUT/'study.json').read_text(encoding='utf-8'))
    checks = []
    def check(name, passed, error=None, tolerance=None):
        row = {'name': name, 'passed': bool(passed)}
        if error is not None: row['error'] = mp.nstr(error, 30)
        if tolerance is not None: row['tolerance'] = str(tolerance)
        checks.append(row)
    def nu_for(label):
        return mp.mpc(0,'0.5') if label=='i/2' else mp.mpf(label)
    for row in study['baseline']:
        label = f"nu={row['nu']},q={row['q']}"
        e_hs = mp.mpf(row['hyper_series_abs_error'])
        e_ha = mp.mpf(row['hyper_angle_abs_error'])
        bound = mp.mpf(row['series_tail_bound'])
        check('series tail '+label, e_hs <= bound + mp.mpf('1e-52'), e_hs, bound+mp.mpf('1e-52'))
        check('geometric integral '+label, e_ha < mp.mpf('1e-48'), e_ha, '1e-48')
        check('real correlator '+label, abs(mp.mpf(row['hyper_imag'])) < mp.mpf('1e-48'))
        if row['nu']=='i/2':
            expected=mp.atanh(mp.mpf(row['q']))/mp.mpf(row['q'])
            err=abs(mp.mpf(row['hyper_real'])-expected)
            check('conformal independent closed form '+label, err < mp.mpf('1e-48'),err,'1e-48')
    for row in study['branches']:
        label=f"q={row['q_real']},t_sign={row['time_sign']},{row['contour']},eps={row['epsilon']}"
        err=mp.mpf(row['elementary_abs_error'])
        check('regulated elementary continuation '+label,err<mp.mpf('1e-48'),err,'1e-48')
        check('time-to-q side '+label,mp.sign(mp.mpf(row['q_imag_regulated']))==row['q_cut_side'])
    for q in ['1.05','1.2','2','5']:
        for sign_t in [-1,1]:
            for contour in ['F12_upper_t','F14_lower_t']:
                rows=[r for r in study['branches'] if r['q_real']==q and r['time_sign']==sign_t and r['contour']==contour]
                errors=[mp.mpf(r['boundary_abs_error']) for r in rows]
                ratios=[errors[i+1]/errors[i] for i in range(len(errors)-1)]
                check(f'linear regulator convergence q={q},t_sign={sign_t},{contour}',
                      all(mp.mpf('.009') < ratio < mp.mpf('.011') for ratio in ratios))
    for label in ['0','0.5','1','2','i/2']:
        nu=nu_for(label)
        check('zero recoil '+label,recoil_hyper(nu,mp.mpf(0))==1 and recoil_angle(nu,0)==1)
        q=mp.mpf('1e-8')
        coefficient=(recoil_hyper(nu,q)-1)/(q*q)
        expected=4*(1+nu*nu)/9
        error=abs(coefficient-expected)
        check('leading small-q coefficient '+label,error < mp.mpf('1e-13'),error,'1e-13')
    # The same fixed physical time and heavy mass must approach the no-recoil value.
    deviations=[]
    for m in ['10000','20000','40000']:
        q=q_from_time(m,'100',mp.mpf(2))
        deviations.append(abs(recoil_hyper(1,q)-1))
    check('fixed-time mass limit',all(deviations[i+1]<deviations[i] for i in range(2)))
    check('leading 1/m^2 scaling',all(abs(deviations[i]/deviations[i+1]-4)<mp.mpf('.01') for i in range(2)))
    # Check the same cut from opposite real-time directions without assuming an i0 sign.
    for sign_t in [-1,1]:
        _, upper=time_branch_q('1.2','0.01',sign_t,1,'1e-8')
        _, lower=time_branch_q('1.2','0.01',sign_t,-1,'1e-8')
        err=abs(recoil_hyper(1,upper)-mp.conj(recoil_hyper(1,lower)))
        check('principal-series conjugate branches t_sign='+str(sign_t),err<mp.mpf('1e-55'),err,'1e-55')
    # Conditional diagonal spectral insertion only: this is not a clock-state model.
    for s_strings in [('0.4','-0.7','0.1','0.2'),('-0.5','0.2','0.1','0.2'),('0','0','0','0')]:
        shifts=list(map(mp.mpf,s_strings))
        s_min=min(shifts)
        for kappa_str in ['-1','0.3','2']:
            kappa=mp.mpf(kappa_str)
            integral=mp.quad(lambda x: mp.exp(-2*mp.pi*x)*mp.exp(-kappa*(x+shifts[0])),
                             [-s_min, -s_min+1, mp.inf])
            formula=mp.exp(2*mp.pi*s_min-kappa*(shifts[0]-s_min))/(2*mp.pi+kappa)
            error=abs(integral-formula)
            check('conditional F4 insertion '+str(s_strings)+',kappa='+kappa_str,
                  error < mp.mpf('1e-55'),error,'1e-55')
            normalized=integral*(2*mp.pi+kappa)
            target=mp.exp(2*mp.pi*s_min-kappa*(shifts[0]-s_min))
            check('conditional F5 normalization '+str(s_strings)+',kappa='+kappa_str,
                  abs(normalized-target)<mp.mpf('1e-55'))
    for name, expected_hash in study['code_sha256'].items():
        actual=hashlib.sha256((ROOT/'code'/name).read_bytes()).hexdigest()
        check('code hash '+name,actual==expected_hash)
    failed=sum(not r['passed'] for r in checks)
    result={'precision_dps':70,'status':'failed' if failed else 'passed',
            'checks':len(checks),'check_details':checks,'passed':sum(r['passed'] for r in checks),
            'failed':failed,'code_sha256':study['code_sha256'],
            'scope':'Formula reproduction, regulator diagnostics, limiting cases, conditional diagonal insertion; no finite-clock dictionary validation.'}
    (OUT/'validation.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(f"{result['passed']} passed; {result['failed']} failed")
    for row in checks:
        if not row['passed']: print(json.dumps(row))
    if result['failed']: raise SystemExit(1)


if __name__=='__main__':
    main()
