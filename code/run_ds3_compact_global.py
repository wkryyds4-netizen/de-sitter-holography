"""New full-kernel samples outside the former small-numerator window."""
import hashlib,json,platform
from pathlib import Path
import mpmath as mp
import sympy as sp
from ds3_compact_global import evaluate,hilbert_tail
from ds3_shadow_duality import physical
from run_ds3_collision import clean


def main():
    root=Path(__file__).resolve().parents[1];base=root/'results/stage4/ds3-compact-global'
    base.mkdir(parents=True,exist_ok=True);raw=base/'runs.jsonl';dest=base/'study.json'
    if raw.exists() or dest.exists():raise RuntimeError('Preserve old output')
    mp.mp.dps=90;rows=[];references=[];products=[]
    families=[('1.4','-.8','1.2','1.7'),('-.9','1.1','2.2','.6')]
    for M in [16,64,256]:
        for family in families:
            for r in ['-.6','0','.6']:
                row=evaluate(M,*family,r,R='1',N=160)
                rows.append(clean(row))
                with raw.open('a',encoding='utf-8') as out:out.write(json.dumps(rows[-1])+'\n');out.flush()
                print('kernel',M,family,r,'bound',row['passed'],'sine err',mp.nstr(row['sine_error'],8),flush=True)
    for index in [0,1,5,9,16,17]:
        row=rows[index];M,u,v,p,q,r=map(mp.mpf,[row[k] for k in ['M','u','v','p','q','r']])
        with mp.workdps(120):
            ref,tail,_=physical(str(M+u),str(M+v+r),list(map(str,[M+u+v,q,M,p])),refinement=4)
            old=mp.mpc(*map(mp.mpf,row['value']));diff=abs(old-ref);allow=mp.mpf('1e-65')
            references.append(clean(dict(index=index,value=ref,tail=tail,difference=diff,rounding_allowance=allow,
                                         passed=diff<=tail+mp.mpf(row['tail'])+allow)))
        print('reference',index,references[-1]['passed'],flush=True)
    for alpha in map(mp.mpf,['0','.5','1','3','8']):
        for n in [1,2,10,100]:
            actual=abs(mp.rf(mp.mpf('.5')+1j*alpha,n)/mp.factorial(n))
            bound=mp.sqrt(mp.cosh(mp.pi*alpha)/(n+1))
            products.append(clean(dict(alpha=alpha,n=n,actual=actual,bound=bound,passed=actual<=bound)))
    paths=['code/ds3_compact_global.py','code/run_ds3_compact_global.py','code/ds3_shadow_duality.py',
           'code/ds3_full_tail_bound.py','code/ds3_transverse_phase.py']
    result=dict(status='passed' if all(x['passed'] for x in rows+references+products) else 'failed',rows=rows,
                references=references,product_checks=products,tails=[clean(hilbert_tail(K,'.2')) for K in [0,2,4,8,12,16]],
                dps=90,reference_dps=120,dependencies=dict(python=platform.python_version(),mpmath=mp.__version__,sympy=sp.__version__),
                scope='Full scalar kernels at 18 finite points and six unmerged references; not full physical quadrature or interval certification. Global limit is an analytic theorem in the note.',
                hashes={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in paths})
    dest.write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(result['status'],flush=True)
    if result['status']!='passed':raise SystemExit(1)


if __name__=='__main__':main()
