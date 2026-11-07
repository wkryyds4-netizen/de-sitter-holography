"""All endpoint orbits of the cubic residual versus either force orientation."""
from pathlib import Path
from itertools import combinations,permutations
from math import factorial
import json,hashlib,platform
import sympy as s
from quartic_orbits import cross,polynomial_dict,accumulate,q
from memory_oriented_marks import graph_oriented

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-mixed-spectrum'


def certificate(reverse):
    ma=list(combinations(range(3),2));mb=[(0,1)]
    D=(1+q*q)*(1+q+q*q)
    wa=[-q*(1+q+q*q),D,-q*(1+q+q*q),-2*(1-q)*(1+q**3)]
    wb=[1+q,-2*(1-q)]
    prod=[[polynomial_dict(a*b) for b in wb] for a in wa]
    records=[];fails=[];count=checks=0
    for k in range(3):
        n=5-2*k
        for aset in combinations(range(n),3-k):
            colors=''.join('A' if j in aset else 'B' for j in range(n))
            for split in range(n+1):
                A=list(range(3));B=list(range(3+split,5+split))
                ext=list(range(3,3+split))+list(range(5+split,5+n))
                ae=[v for v,c in zip(ext,colors) if c=='A'];be=[v for v,c in zip(ext,colors) if c=='B']
                seed=list(zip(A[:k],B[:k]))+list(zip(A[k:],ae))+list(zip(B[k:],be))
                words=set()
                for ap in permutations(A):
                    for bp in permutations(B):
                        pm=dict(zip(A,ap))|dict(zip(B,bp))
                        words.add(tuple(sorted(tuple(sorted((pm.get(a,a),pm.get(b,b)))) for a,b in seed)))
                assert len(words)==12//factorial(k)
                linear={};quadratic={};count+=len(words)
                for pa in sorted(words):
                    pi={a:b for a,b in pa}|{b:a for a,b in pa}
                    cr=sum(cross(a,b,c,d) for j,(a,b) in enumerate(pa) for c,d in pa[j+1:])
                    al=[2-4*cross(A[a],pi[A[a]],A[b],pi[A[b]]) for a,b in ma]+[1]
                    be=[2-4*cross(B[0],pi[B[0]],B[1],pi[B[1]]),1]
                    for i in range(4):
                        for j in range(2):accumulate(quadratic,prod[i][j],al[i]*be[j],cr)
                    for i,(a,b) in enumerate(ma):
                        dd=(-4*int(pi[A[a]]==B[0] and pi[A[b]]==B[1]) if reverse
                            else 4*int(pi[A[a]]==B[1] and pi[A[b]]==B[0]))
                        g=graph_oriented(5+n,pa,((A[a],A[b]),(B[0],B[1])),(False,reverse))
                        assert (g['crossings'],g['alpha'],g['beta'],g['delta'])==(cr,al[i],be[0],dd)
                        checks+=1;accumulate(linear,prod[i][0],dd,cr)
                linear={j:v for j,v in linear.items() if v};quadratic={j:v for j,v in quadratic.items() if v}
                row=dict(bridges=k,colors=colors,split=split,words=len(words),linear=linear,quadratic=quadratic)
                records.append(row)
                if linear or quadratic:fails.append(row)
    return dict(reverse=reverse,templates=len(records),words=count,minor_checks=checks,
                nonzero_templates=fails,records=records)


def main():
    OUT.mkdir(exist_ok=True);target=OUT/'cubic-force-cross.json'
    if target.exists():raise FileExistsError(target)
    rows=[certificate(False),certificate(True)]
    names=['code/memory_mixed_cross.py','code/quartic_orbits.py','code/memory_oriented_marks.py']
    data=dict(status='passed' if all(not r['nonzero_templates'] for r in rows) else 'nonzero cross',
              certificates=rows,environment=dict(python=platform.python_version(),sympy=s.__version__),
              input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names})
    target.write_text(json.dumps(data,indent=2)+'\n')
    print([(r['reverse'],r['templates'],r['words'],len(r['nonzero_templates'])) for r in rows],flush=True)


if __name__=='__main__':main()
