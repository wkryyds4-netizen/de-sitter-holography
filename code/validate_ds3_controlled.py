"""Independent hyper/parameter-derivative comparisons for analytic tail bounds."""
import hashlib
import json
from pathlib import Path
import mpmath as mp
import sympy as sp
from ds3_controlled_omega import controlled
from ds3_racah_collision import omega_collision_derivative
from global_racah import omega_hyper
from run_ds3_collision import clean


def main():
    root=Path(__file__).resolve().parents[1];dest=root/'results/stage4/ds3-collision/controlled-validation.json'
    if dest.exists():raise RuntimeError('Preserve prior run')
    mp.mp.dps=100;rows=[]
    for M in [16,32,64]:
        c=mp.mpf(M);a=c+mp.mpf('.1');b=c+mp.mpf('.3');d0=c-mp.mpf('.2');p=mp.mpf('.6')
        hs=[(1+1j*x)/2 for x in [a,p,c,p]];H=(1-1j*b)/2
        for lam in [0,1]:
            r=mp.mpf(lam)/mp.log(M) if lam else mp.mpf(0);hp=(1+1j*(d0+r))/2
            reference=omega_collision_derivative(H,hs) if not lam else omega_hyper(H,hp,hs)[0]
            for N in [32,128]:
                val,detail=controlled(M,'.3','-.2','.6','.6',r,N)
                err=abs(val-reference);tail=detail['series_tail_bound']
                rows.append(dict(M=M,lambda_value=lam,N=N,value=val,reference=reference,error=err,
                                 tail=tail,tail_passed=err<=tail+mp.mpf('1e-85'),
                                 leading_error=abs(reference-detail['leading']),leading_bound=detail['leading_bound'],
                                 leading_passed=abs(reference-detail['leading'])<=detail['leading_bound']))
                print(M,lam,N,mp.nstr(err,5),mp.nstr(tail,5),flush=True)
    x=sp.symbols('x',positive=True)
    algebra=sp.expand(x**3-(x*x-x+sp.Rational(1,2))*(x+1))
    assert algebra==(x-1)/2
    paths=['code/validate_ds3_controlled.py','code/ds3_controlled_omega.py','code/ds3_collision_asymptotic.py',
           'code/ds3_racah_collision.py','code/global_racah.py']
    data=dict(status='passed' if all(r['tail_passed'] and r['leading_passed'] for r in rows) else 'failed',
              dps=100,rounding_allowance='1e-85 absolute, empirical rather than interval',rows=rows,
              ratio_symbolic_difference=str(algebra),
              hashes={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in paths})
    dest.write_text(json.dumps(clean(data),indent=2),encoding='utf-8')


if __name__=='__main__':main()
