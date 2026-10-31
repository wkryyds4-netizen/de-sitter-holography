"""Small nondegenerate point validation for the exact observer crossing kernel.

No multidimensional integration, removable-pole extrapolation, or interval
certification. Failed diagnostics are saved and produce a nonzero exit code.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import sys
import traceback

import mpmath
from observer_crossing import CrossingEvaluator

ROOT = Path(__file__).resolve().parents[1]
POINTS = [
    {'nu':'0.6','labels':['0.35','0.7','1.1','1.45']},
    {'nu':'1.0','labels':['0.4','1.25','0.8','1.7']},
]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT/'results/stage3/fourpoint/crossing-validation.json')
    args = parser.parse_args()
    low_dps,high_dps,guard = 25,40,10
    low = CrossingEvaluator(low_dps,guard)
    high = CrossingEvaluator(high_dps,guard)
    c = high.ctx
    tolerance = c.mpf('1e-20')
    checks = []
    records = []
    failures = []

    def check(name,residual,limit=tolerance):
        residual = c.mpf(residual)
        checks.append({'name':name,'passed':bool(c.isfinite(residual) and residual <= limit),
                       'residual':c.nstr(residual,25),'limit':c.nstr(limit,10)})

    def relative(x,y):
        x,y = c.mpc(x),c.mpc(y)
        scale = max(abs(x),abs(y),c.mpf('1e-100'))
        return abs(x-y)/scale

    def complex_record(x,digits=50):
        x = c.mpc(x)
        return {'real':c.nstr(c.re(x),digits),'imag':c.nstr(c.im(x),digits)}

    def record(result):
        return {'labels':[str(v) for v in result.labels],
                'nu':str(result.matter_nu),'requested_dps':result.requested_dps,
                'working_dps':result.working_dps,'elapsed_seconds':result.elapsed_seconds,
                'raw_s':complex_record(result.raw_s),
                'gamma_product':complex_record(result.gamma_product),
                'symmetric_s':complex_record(result.symmetric_s),
                'hh_cross':complex_record(result.hh_cross),
                'g_terms':[complex_record(v) for v in result.g_terms],
                'cancellation_ratio':c.nstr(result.cancellation_ratio,25),
                'hyper_calls':result.hyper_calls,
                'hyper_excess_errors':[c.nstr(v,10) for v in result.hyper_excess_errors]}

    for index,point in enumerate(POINTS):
        nu,labels = point['nu'],point['labels']
        name = f'point_{index+1}'
        print(f'{name}: evaluating D4 orbit and precision refinement',flush=True)
        point_record = {'input':point,'orbit':[]}
        try:
            first = low.evaluate(nu,*labels)
            refined = high.evaluate(nu,*labels)
            point_record['base_low'] = record(first)
            point_record['base_high'] = record(refined)
            check(name+'/precision_symmetric_s',relative(first.symmetric_s,refined.symmetric_s))
            check(name+'/precision_hh_cross',relative(first.hh_cross,refined.hh_cross))
            check(name+'/reality_high',abs(c.im(refined.symmetric_s))/max(abs(refined.symmetric_s),c.mpf('1e-100')))
            direct = high.gamma_product_expanded(refined.labels,refined.matter_nu)
            check(name+'/A14_direct_gamma_product',relative(direct,refined.gamma_product),c.mpf('1e-40'))
            all_results = {}
            for reflected in (False,True):
                orientation = list(reversed(labels)) if reflected else labels
                for shift in range(4):
                    perm = orientation[shift:]+orientation[:shift]
                    key = tuple(perm)
                    result = first if key == tuple(labels) else low.evaluate(nu,*perm)
                    all_results[key] = result
                    tag = f'{name}/reflection_{int(reflected)}/rotation_{shift}'
                    check(tag+'/symmetric_s',relative(result.symmetric_s,first.symmetric_s))
                    check(tag+'/reality',abs(c.im(result.symmetric_s))/max(abs(result.symmetric_s),c.mpf('1e-100')))
                    check(tag+'/hyper_excess',max(result.hyper_excess_errors),c.mpf('1e-30'))
                    # Complete X is not cyclic: its explicit odd-position cosh
                    # weights must compensate the symmetric part.
                    a,b,z,d = [c.mpf(v) for v in labels]
                    p,q,r,s = [c.mpf(v) for v in perm]
                    factor = c.cosh(c.pi*a)*c.cosh(c.pi*z)/(c.cosh(c.pi*p)*c.cosh(c.pi*r))
                    check(tag+'/full_cross_weighted_symmetry',relative(first.hh_cross,factor*result.hh_cross))
                    point_record['orbit'].append({'reflected':reflected,'shift':shift,
                                                  'permutation':perm,'result':record(result)})
            reversed_result = all_results[tuple(reversed(labels))]
            a,b,z,d = [c.mpf(v) for v in labels]
            ratio = c.cosh(c.pi*a)*c.cosh(c.pi*z)/(c.cosh(c.pi*d)*c.cosh(c.pi*b))
            check(name+'/source_6_27',relative(first.hh_cross,ratio*reversed_result.hh_cross))
            point_record['source_6_27_factor'] = c.nstr(ratio,40)
            point_record['full_cross_plain_reversal_residual'] = c.nstr(relative(first.hh_cross,reversed_result.hh_cross),25)
        except Exception as exc:
            failure = {'point':point,'type':type(exc).__name__,'message':str(exc),
                       'traceback':traceback.format_exc()}
            failures.append(failure)
            point_record['failure'] = failure
        records.append(point_record)

    for label,args_bad in [('opposite_labels',('0.6','0.4','0.7','0.4','1.1')),
                           ('zero_mass_label',('0.6','0','0.7','0.9','1.1'))]:
        try:
            low.evaluate(*args_bad)
        except ValueError:
            checks.append({'name':'domain/'+label,'passed':True,'expected':'ValueError'})
        except Exception as exc:
            checks.append({'name':'domain/'+label,'passed':False,'actual':repr(exc)})
        else:
            checks.append({'name':'domain/'+label,'passed':False,'actual':'accepted'})

    passed = sum(row['passed'] for row in checks)
    report = {
        'schema':'observer-crossing-validation-v1',
        'created_utc':datetime.now(timezone.utc).isoformat(),
        'model':{'spacetime':'dS2','L':1,'G_N':'strict zero limit',
                 'matter':'two independent principal-series species with equal nu',
                 'kernel':'unnormalized HH crossing density in da db dc dd',
                 'external_mass':'d','source':'https://arxiv.org/html/2411.08090v1',
                 'source_equations':['6.26','6.27','A.14','B.21','B.22','B.23','B.24']},
        'numerics':{'library':'mpmath','version':mpmath.__version__,
                    'python':sys.version,'platform':platform.platform(),
                    'low_requested_dps':low_dps,'high_requested_dps':high_dps,
                    'guard_digits':guard,'sum_method':'r+s+e','strict_convergence':True,
                    'tolerance':str(tolerance),'interval_certification':False},
        'provenance':{str(p.relative_to(ROOT)):digest(p) for p in [
            ROOT/'code/observer_crossing.py',Path(__file__).resolve(),
            ROOT/'sources/kolchmeyer-liu-2411.08090v1.html',
            ROOT/'literature/stage3-fourpoint-source-map.md']},
        'summary':{'checks':len(checks),'passed':passed,'failed':len(checks)-passed,
                   'evaluation_exceptions':len(failures),
                   'all_passed':passed==len(checks) and not failures},
        'checks':checks,'points':records,'exceptions':failures,
        'limitations':['Nondegenerate positive labels only; coincident-label limits are unsupported.',
                       'Only two D4 orbits; precision refinement and symmetry checks are not interval bounds.',
                       'The common normalization is taken from the cited source, not fixed independently by symmetry.',
                       'No four-mass wavepacket integration or proof of its absolute convergence.',
                       'All raw complex values are retained; no forced real-part fallback.']}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report['summary']),flush=True)
    print(args.output,flush=True)
    return 0 if report['summary']['all_passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
