"""Supplement sparse eight-slot random words with signed nonzero CAR words."""
from pathlib import Path
from fractions import Fraction
from itertools import combinations
import json,hashlib
import numpy as np
from complex_syk_dictionary import fermions,creators
from check_singlet_spectral_kernel import fixed_word

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-singlet-kernel/nonzero-words.json'
SOURCE='results/susskind-singlet-double-scale/word-checks.json'

def main():
    if OUT.exists():raise FileExistsError(OUT)
    old=json.loads((ROOT/SOURCE).read_text());records=[];nonzero=negative=0
    cached={}
    for row in old['records']:
        if Fraction(row['exact_values'][3])==0:continue
        n,p=row['N'],row['p'];d=1<<n
        if (n,p) not in cached:
            cs=fermions(n);ss=list(combinations(range(n),p))
            cached[n,p]=(cs,ss,dict(zip(ss,creators(cs,p))),np.diag([s.bit_count() for s in range(d)]))
        cs,subsets,bs,Q=cached[n,p]
        pairing=tuple(map(tuple,row['pairing']))+((6,7),)
        marks=tuple(map(tuple,row['marks']));pi={a:b for a,b in pairing}|{b:a for a,b in pairing}
        for extra in subsets:
            sets=tuple(map(tuple,row['creator_sets']))+(extra,extra)
            ops=[bs[sets[a]]@bs[sets[pi[a]]].conj().T for a in range(8)]
            def deriv(x,y):
                z=x@y
                return sum(c.conj().T@z@c for c in cs)-(Q-2*p*np.eye(d))@z
            values=[]
            for enabled in [(),(0,),(1,),(0,1)]:
                starts={marks[j][0] for j in enabled};z=np.eye(d,dtype=complex);i=0
                while i<8:
                    if i in starts:z=z@deriv(ops[i],ops[i+1]);i+=2
                    else:z=z@ops[i];i+=1
                v=np.trace(z)
                assert v.imag==0 and v.real==round(v.real)
                values.append(Fraction(int(v.real),d))
            assert tuple(values)==fixed_word(n,p,sets,pairing,marks)
            nonzero+=int(values[3]!=0);negative+=int(values[3]<0)
            records.append(dict(N=n,p=p,sets=sets,pairing=pairing,marks=marks,values=list(map(str,values))))
    assert nonzero>=50 and negative>0
    data=dict(status='passed',selection='Extend every prior nonzero six-slot CAR word by an adjacent paired diagonal p-body projector; all p-subsets used.',
              word_cases=len(records),scalar_comparisons=4*len(records),nonzero=nonzero,negative=negative,records=records,
              input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in
              [SOURCE,'code/check_singlet_kernel_nonzero.py','code/check_singlet_spectral_kernel.py']})
    with OUT.open('x') as f:json.dump(data,f,indent=2);f.write('\n')
    print(json.dumps({k:v for k,v in data.items() if k not in ['records','input_sha256']}))

if __name__=='__main__':main()
