"""Independent exact finite CAR check of the multi-mark determinant formula."""
from pathlib import Path
from itertools import combinations
from fractions import Fraction
import hashlib,json,platform
import numpy as np
from singlet_coefficient import determinant,pairings
from singlet_multimark import nmul
from complex_syk_dictionary import fermions,creators

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-singlet-multimark'


def local_coeff(n,p,sets,pairing,marks):
    length=len(sets);pi={a:b for a,b in pairing}|{b:a for a,b in pairing}
    product={0:1}
    for mode in range(n):
        slots=[i for i in range(length) if mode in sets[i]];vals=[]
        for mask in range(1<<len(marks)):
            vals.append(determinant([[(1 if a<=pi[b] else -1)-2*sum(int((pi[b],a)==m) for j,m in enumerate(marks) if mask>>j&1) for b in slots] for a in slots]))
        for j in range(len(marks)):
            for mask in range(1<<len(marks)):
                if mask>>j&1:vals[mask]-=vals[mask^(1<<j)]
        product=nmul(product,{k:v for k,v in enumerate(vals) if v})
    return [Fraction((-1)**(length//2*p)*product.get(mask,0),2**(length*p)) for mask in range(1<<len(marks))]


def main():
    target=OUT/'finite-car.json'
    if target.exists():raise FileExistsError(target)
    rng=np.random.default_rng(110609);rows=[];nonzero=negative=0
    for n,p,m in [(3,1,3),(4,2,3),(3,1,4),(4,2,4),(3,1,6)]:
        cs=fermions(n);ss=list(combinations(range(n),p));bs=dict(zip(ss,creators(cs,p)))
        d=1<<n;Q=np.diag([x.bit_count() for x in range(d)])
        marks=[(2*i,2*i+1) for i in range(m)]
        allowed=[r for r in pairings(tuple(range(2*m))) if not any(x in r for x in marks)]
        for case in range(32):
            pairing=allowed[int(rng.integers(len(allowed)))];pi={a:b for a,b in pairing}|{b:a for a,b in pairing}
            # Include coherent repeated mode sets, not only sparse random zeros.
            sets=[ss[int(rng.integers(len(ss)))]]*(2*m) if case<8 else [ss[int(rng.integers(len(ss)))] for _ in range(2*m)]
            ops=[bs[sets[a]]@bs[sets[pi[a]]].conj().T for a in range(2*m)]
            plain=[ops[a]@ops[b] for a,b in marks]
            marked=[sum(c.conj().T@v@c for c in cs)-(Q-2*p*np.eye(d))@v for v in plain]
            actual=[]
            for mask in range(1<<m):
                mat=np.eye(d,dtype=complex)
                for j in range(m):mat=mat@(marked[j] if mask>>j&1 else plain[j])
                value=np.trace(mat)
                assert value.imag==0 and value.real==round(value.real)
                actual.append(Fraction(int(value.real),d))
            expected=local_coeff(n,p,sets,pairing,marks)
            assert actual==expected,(n,p,m,sets,pairing,actual,expected)
            nonzero+=int(actual[-1]!=0);negative+=int(actual[-1]<0)
            rows.append(dict(N=n,p=p,marks=marks,sets=sets,pairing=pairing,values=list(map(str,actual))))
        print(n,p,m,'passed',flush=True)
    assert nonzero>=40
    inputs=['code/check_multimark_car.py','code/singlet_multimark.py','code/singlet_coefficient.py','code/complex_syk_dictionary.py']
    data=dict(status='passed',seed=110609,words=len(rows),comparisons=sum(len(r['values']) for r in rows),full_mark_nonzero=nonzero,full_mark_negative=negative,
        scope='finite fixed-word CAR versus determinant; checks signs and all mark subsets, not finite-size convergence',
        environment=dict(python=platform.python_version(),numpy=np.__version__),records=rows,
        input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs})
    target.write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps({k:data[k] for k in ['status','words','comparisons','full_mark_nonzero','full_mark_negative']}))

if __name__=='__main__':main()
