"""Exact finite combinatorics with independent 80/120-digit coefficient evaluation."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,platform
import mpmath as mp
from memory_two_charge_filter import filters,choose
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-two-charge'

def main():
    target=OUT/'limit-checks.json'
    if target.exists():raise FileExistsError(target)
    records=[]
    for lam in [F(1,10),F(1,4),F(1,2)]:
        for p in [16,32,64,128]:
            N=int(p*p/lam);assert F(p*p,N)==lam and N%2==0
            m=N//2;C=choose(m,p)*choose(N-m+p,p)
            for u,v,S in [(0,0,0),(2,1,3)]:
                gp,gm=filters(N,p,m,p-u,p-v,2*p-u-v-S)['g']
                exactA=gm/C;exactH=(gp+(m-2*p)*gm)/C
                values=[]
                for precision in [80,120]:
                    mp.mp.dps=precision
                    cv=lambda f:mp.mpf(f.numerator)/f.denominator
                    l=cv(lam);q=mp.exp(-4*l)
                    gA=cv(exactA);gH=cv(exactH)
                    Alead=(3-4*q+q*q)/2;Hlead=2*p*(1-q)**2
                    Hfinite=l*(8*S*q*q-32*l*q*q+16*l*q+8*(q*q-q)*(u+v)+9*q*q-10*q+1)
                    values.append((gA,gH,Hfinite,gH-Hlead))
                for low,high in zip(values[0],values[1]):assert abs(low-high)<mp.mpf('1e-70')*max(1,abs(high))
                k=choose(m-(p-u),p)*choose(N-m+p-(p-u),p)
                actual_k=cv(F(k,C));approx_k=q*(1+(4*l*u-4*l*l)/p)
                records.append(dict(N=N,p=p,m=m,lambda_value=str(lam),input_deficits=[u,v],output_deficit=S,
                    gA_over_C=mp.nstr(gA,60),gA_limit=mp.nstr(Alead,60),
                    gA_absolute_deviation=mp.nstr(gA-Alead,40),
                    gH_over_C_minus_leading=mp.nstr(gH-Hlead,60),gH_finite_limit=mp.nstr(Hfinite,60),
                    gH_absolute_deviation=mp.nstr(gH-Hlead-Hfinite,40),
                    kappa_first_order_remainder_times_p_squared=mp.nstr(p*p*(actual_k-approx_k),40)))
    names=['code/check_memory_filter_limit.py','code/memory_two_charge_filter.py',
        'results/susskind-memory-two-charge/symbolic-filter.json']
    result=dict(status='passed',records=records,precisions=[80,120],
        scope='Finite deviations retained. Decimal precision agreement is not an interval enclosure or an operator-limit proof.',
        environment=dict(python=platform.python_version(),mpmath=mp.__version__),
        input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names})
    target.write_text(json.dumps(result,indent=2)+'\n')
    print('passed',len(records),'finite cases; last',records[-1])

if __name__=='__main__':main()
