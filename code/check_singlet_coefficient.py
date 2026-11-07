"""Compare marked local determinants with direct CAR matrices, before averaging."""
from fractions import Fraction
from itertools import combinations
import numpy as np
from complex_syk_dictionary import fermions,creators
from singlet_coefficient import fixed_word,pairings,ALL_MARKS


def check_words():
    rng=np.random.default_rng(110506)
    records=[]
    nonzero=negative=0
    for n,p in [(4,1),(5,2),(6,3)]:
        cs=fermions(n)
        subsets=list(combinations(range(n),p))
        creators_by_set=dict(zip(subsets,creators(cs,p)))
        d=1<<n
        q=np.diag([s.bit_count() for s in range(d)])
        def derivative(x,y):
            xy=x@y
            return sum(c.conj().T@xy@c for c in cs)-(q-2*p*np.eye(d))@xy
        for pairing in pairings(tuple(range(6))):
            pi={a:b for a,b in pairing}|{b:a for a,b in pairing}
            for draw in range(4):
                sets=[subsets[i] for i in rng.integers(0,len(subsets),size=6)]
                ops=[creators_by_set[sets[a]]@creators_by_set[sets[pi[a]]].conj().T for a in range(6)]
                for marks,_ in ALL_MARKS:
                    direct=[]
                    for enabled in [(),(0,),(1,),(0,1)]:
                        starts={marks[j][0] for j in enabled}
                        result=np.eye(d,dtype=complex)
                        idx=0
                        while idx<6:
                            if idx in starts:
                                result=result@derivative(ops[idx],ops[idx+1]);idx+=2
                            else:
                                result=result@ops[idx];idx+=1
                        value=np.trace(result)
                        assert value.imag==0 and value.real==round(value.real)
                        direct.append(Fraction(int(value.real),d))
                    reduced=fixed_word(n,p,sets,pairing,marks)
                    assert tuple(direct)==reduced,(n,p,sets,pairing,marks,direct,reduced)
                    nonzero+=int(direct[3]!=0);negative+=int(direct[3]<0)
                    records.append(dict(N=n,p=p,pairing=pairing,creator_sets=sets,marks=marks,
                                        exact_values=[str(v) for v in direct]))
    assert nonzero>20 and negative>0
    return dict(status='passed',seed=110506,word_cases=len(records),scalar_comparisons=4*len(records),
                nonzero_mixed_words=nonzero,negative_mixed_words=negative,records=records)
