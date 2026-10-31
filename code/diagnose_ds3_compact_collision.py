"""Preserve and diagnose the failed unmerged reference; exact rational labels."""
import hashlib,json,platform
from pathlib import Path
import mpmath as mp
import sympy as sp
from ds3_shadow_duality import physical
from run_ds3_collision import clean

def main():
    root=Path(__file__).resolve().parents[1]
    base=root/'results/stage4/ds3-compact-global'
    dest=base/'collision-diagnosis.json'
    if dest.exists():raise RuntimeError('Preserve existing diagnosis')
    src=base/'study.json';study=json.loads(src.read_text());row=study['rows'][16]
    # Reproduce the exact driver ordering, including 90-digit input parsing.
    mp.mp.dps=90
    M,u,v,p,q,r=[mp.mpf(row[k]) for k in ['M','u','v','p','q','r']]
    with mp.workdps(120):
        rounded=list(map(str,[M+u,M+v+r,M+u+v,q,M,p]))
    b,d,a,qq,c,pp=map(sp.Rational,rounded)
    residual=sp.factor(a+c-b-d)
    exact=['2551/10','2571/10','1281/5','3/5','256','11/5']
    assert sp.Rational(exact[2])+sp.Rational(exact[4])-sp.Rational(exact[0])-sp.Rational(exact[1])==0
    results=[]
    for mode,labels,dps in [('rounded',rounded,120),('rounded',rounded,210),('exact',exact,120),('exact',exact,170)]:
        with mp.workdps(dps):
            val,tail,details=physical(labels[0],labels[1],labels[2:],refinement=4)
            old=mp.mpc(*map(mp.mpf,row['value']))
            err=abs(val-old);allow=mp.mpf('1e-65')
            y=details['omega'][1]['details']
            result=clean(dict(mode=mode,dps=dps,value=val,tail=tail,difference=err,
                relative_difference=err/max(abs(old),mp.mpf('1e-1000')),
                collision=y['collision'],z=y.get('z'),pair_condition=y.get('pair_condition'),
                passed=err<=tail+mp.mpf(row['tail'])+allow))
            results.append(result)
            print(mode,dps,'collision',y['collision'],'difference',mp.nstr(err,8),'passed',result['passed'],flush=True)
    paths=['code/diagnose_ds3_compact_collision.py','code/ds3_shadow_duality.py','code/ds3_full_tail_bound.py',str(src.relative_to(root))]
    report=dict(status='passed' if not results[0]['passed'] and all(x['passed'] for x in results[1:]) else 'failed',
        index=16,rounded_labels=rounded,exact_labels=exact,exact_collision_residual=str(residual),
        rows=results,rounding_allowance='1e-65',dependencies=dict(python=platform.python_version(),mpmath=mp.__version__,sympy=sp.__version__),
        scope='Reproduced failed reference plus precision and exact-label controls. No interval rounding certification.',
        hashes={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in paths})
    dest.write_text(json.dumps(report,indent=2),encoding='utf-8')
    if report['status']!='passed':raise SystemExit(1)

if __name__=='__main__':main()
