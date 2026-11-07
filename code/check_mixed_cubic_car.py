"""Finite mixed quadratic/cubic CAR versus local determinant, exact dyadic values."""
from pathlib import Path
from itertools import combinations,product
from fractions import Fraction
import json,hashlib,platform
import numpy as np
from complex_syk_dictionary import fermions,creators
from singlet_coefficient import pairings
from check_multimark_car import local_coeff
from mixed_cubic_vertices import blocks
ROOT=Path(__file__).resolve().parent.parent

def main():
    target=ROOT/'results/susskind-mixed-cubic/finite-car.json'
    if target.exists():raise FileExistsError(target)
    rng=np.random.default_rng(110619);rows=[];checks=nonzero=negative=0
    for n,p in [(3,1),(4,2)]:
        cs=fermions(n);ss=list(combinations(range(n),p));bs=dict(zip(ss,creators(cs,p)))
        dim=1<<n;Q=np.diag([i.bit_count() for i in range(dim)])
        def D(v):return Q@v-sum(c.conj().T@v@c for c in cs)
        def block(v):
            a,b,*rest=v;ab=a@b;k12=2*p*ab-D(ab)
            if not rest:return [ab,k12]
            c=rest[0];abc=ab@c;k23=a@(2*p*b@c-D(b@c));k12=k12@c
            return [abc,k12,k23,3*p*abc-D(abc)-k12-k23]
        for word in ['XRR','XRXR','XXRR','RRRR']:
            length,vs,owners=blocks(word)
            allowed=[pa for pa in pairings(tuple(range(length))) if not any(owners[a]==owners[b] for a,b in pa)]
            for case in range(8):
                pa=allowed[int(rng.integers(len(allowed)))];pi={a:b for a,b in pa}|{b:a for a,b in pa}
                sets=[ss[int(rng.integers(len(ss)))]]*length if case<3 else [ss[int(rng.integers(len(ss)))] for _ in range(length)]
                ops=[bs[sets[a]]@bs[sets[pi[a]]].conj().T for a in range(length)]
                matrices=[block([ops[a] for a in v['slots']]) for v in vs]
                choices=[[None]+v['marks'] for v in vs];values=[]
                for ids in product(*(range(len(c)) for c in choices)):
                    mat=np.eye(dim,dtype=complex);marks=[]
                    for j,k in enumerate(ids):
                        mat=mat@matrices[j][k]
                        if k:marks.append(choices[j][k])
                    v=np.trace(mat);assert v.imag==0 and v.real==round(v.real)
                    actual=Fraction(int(v.real),dim);expected=local_coeff(n,p,sets,pa,marks)[-1]
                    assert actual==expected,(n,p,word,case,ids,actual,expected)
                    checks+=1;nonzero+=int(actual!=0);negative+=int(actual<0)
                    values.append(str(actual))
                rows.append(dict(N=n,p=p,word=word,sets=sets,pairing=pa,values=values))
            print(n,p,word,'passed',flush=True)
    assert nonzero>100 and negative>0
    inputs=['code/check_mixed_cubic_car.py','code/check_multimark_car.py','code/mixed_cubic_vertices.py','code/complex_syk_dictionary.py','code/singlet_coefficient.py','code/singlet_multimark.py']
    data=dict(status='passed',seed=110619,words=len(rows),comparisons=checks,nonzero=nonzero,negative=negative,
        scope='Exact finite word and mark-subset CAR identities, not finite-size asymptotics',records=rows,
        environment=dict(python=platform.python_version(),numpy=np.__version__),
        input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs})
    target.write_text(json.dumps(data,indent=2)+'\n')
    print({k:data[k] for k in ['status','words','comparisons','nonzero','negative']})
if __name__=='__main__':main()
