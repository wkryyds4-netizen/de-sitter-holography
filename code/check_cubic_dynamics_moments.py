"""All three Wick pairings for both two-energy degree-two moments."""
from fractions import Fraction as F
from pathlib import Path
import json,hashlib
import numpy as np
import mpmath as mp
from singlet_wick_moments import gue_basis
from cubic_linear_channels import coeff
from cubic_linear_dynamics import finite,weights,b,ds
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-cubic-dynamics'
def main():
    path=OUT/'two-energy-validation.json'
    if path.exists():raise FileExistsError(path)
    rows=[];maxtrace=0
    for N,p,m in [(3,1,2),(4,2,2),(5,2,3)]:
        cs,vs,rs=gue_basis(N,p);dim=1<<N;charge=np.array([a.bit_count() for a in range(dim)]);Q=np.diag(charge)
        D=lambda Y:Q@Y-sum(c.conj().T@Y@c for c in cs)
        def t2(H):
            A=2*p*H@H-D(H@H)
            return 6*p*H@H@H-2*D(H@H@H)-3*(H@A+A@H)
        G=[sum(int(r)*(t2(V+W)+t2(V-W)-2*t2(V)) for W,r in zip(vs,rs)) for V in vs]
        avg=sum(w*coeff(N,p,m,a) for a,w in enumerate(weights(N,m,p)))
        ix=np.flatnonzero(charge==m);h=vs[:,ix[:,None],ix]
        k=np.array([avg.denominator*g[np.ix_(ix,ix)]-8*avg.numerator*v for g,v in zip(G,h)])
        scale=8*avg.denominator;s20=[0,0,0];s11=[0,0,0];norm=0
        def trace_integer(A):
            nonlocal maxtrace
            z=np.trace(A);maxtrace=max(maxtrace,abs(z));assert maxtrace<2**40
            assert z.real==round(z.real) and z.imag==round(z.imag)
            return complex(int(z.real),int(z.imag))
        for a in range(len(h)):
            norm+=int(rs[a])*trace_integer(k[a]@k[a])
            for c in range(len(h)):
                mul=int(rs[a]*rs[c])
                words20=[k[a]@k[a]@h[c]@h[c],k[a]@k[c]@h[a]@h[c],k[a]@k[c]@h[c]@h[a]]
                words11=[k[a]@h[c]@k[a]@h[c],k[a]@h[a]@k[c]@h[c],k[a]@h[c]@k[c]@h[a]]
                for j in range(3):
                    s20[j]+=mul*trace_integer(words20[j]);s11[j]+=mul*trace_integer(words11[j])
        assert all(z.imag==0 for z in s20+s11+[norm])
        assert max(abs(z) for z in s20+s11+[norm])<2**53
        norm=F(int(norm.real),2*scale**2*b(N,m))
        v20=[F(int(z.real),4*scale**2*b(N,m)) for z in s20]
        v11=[F(int(z.real),4*scale**2*b(N,m)) for z in s11]
        formula=finite(N,p,m)
        assert norm==formula['norm'] and sum(v20)==formula['B20'] and sum(v11)==formula['B11']
        rows.append(dict(N=N,p=p,m=m,norm=str(norm),B20_pairings=list(map(str,v20)),B11_pairings=list(map(str,v11))))
    aggregate=[]
    for lam in ['0.1','0.25','0.5']:
        for p in [4,8]:
            N=int(mp.mpf(p*p)/mp.mpf(lam));m=N//2;ex=finite(N,p,m);num=ds(lam,p,100)
            cv=lambda x:mp.mpf(x.numerator)/x.denominator
            exact=[ex['M2']/ex['norm']/ex['C']/F(lam),ex['B20']/ex['C']/ex['norm'],ex['B11']/ex['C']/ex['norm']]
            fields=['normalized_mu2','B20_over_Cnorm','B11_over_Cnorm']
            errors=[abs(cv(x)-mp.mpf(num[f])) for x,f in zip(exact,fields)]
            assert max(errors)<mp.mpf('1e-60')
            aggregate.append(dict(N=N,p=p,lambda_value=lam,max_absolute_error=mp.nstr(max(errors),20)))
    names=['code/check_cubic_dynamics_moments.py','code/cubic_linear_dynamics.py','code/cubic_linear_channels.py','code/cubic_linear_limit.py','code/singlet_wick_moments.py','code/complex_syk_dictionary.py']
    path.write_text(json.dumps(dict(status='passed',all_pairings=rows,aggregate_vs_exact=aggregate,
      max_integer_trace_abs=maxtrace,input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names}),indent=2)+'\n')
    print('passed: 3 direct two-energy sectors, 6 exact aggregate comparisons')
if __name__=='__main__':main()
