"""Exhaustive local clusters for arbitrary H/X/R words, independent weights evaluation."""
from pathlib import Path
from itertools import combinations,permutations
import hashlib,json,platform
import sympy as s
from singlet_coefficient import pairings
from singlet_spectral_graphs import graph
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-mixed-cubic'
q,l=s.symbols('q lambda')
D=(1+q*q)*(1+q+q*q);B=2*(1-q)*(1+q**3)

def choices(kind,slots):
    if kind=='X':return [((slots[0],slots[1]),1+q),(None,-2*l*(1-q))]
    return [((slots[0],slots[1]),-q*(1+q+q*q)),
            ((slots[1],slots[2]),-q*(1+q+q*q)),((slots[0],slots[2]),D),(None,-l*B)]

def value(pa,kinds,blocks):
    length=2*len(pa)
    if len(kinds)==1:
        ans=0
        for mark,c in choices(kinds[0],blocks[0]):
            g=graph(length,pa,(mark,None))
            ans+=c*q**g['crossings']*(1 if mark is None else l*g['alpha'])
        return ans
    ans=0
    for ma,ca in choices(kinds[0],blocks[0]):
        for mb,cb in choices(kinds[1],blocks[1]):
            g=graph(length,pa,(ma,mb))
            if ma is None and mb is None:v=1
            elif ma is None:v=l*g['beta']
            elif mb is None:v=l*g['alpha']
            else:v=l*g['delta']+l*l*g['alpha']*g['beta']
            ans+=ca*cb*q**g['crossings']*v
    return ans

def main():
    target=OUT/'cluster-certificate.json'
    if target.exists():raise FileExistsError(target)
    rows=[]
    for kind,n in [('X',2),('R',3)]:
        block=list(range(n));total=0;count=0
        for perm in permutations(range(n,2*n)):
            pa=tuple(zip(block,perm));total+=value(pa,[kind],[block]);count+=1
        total=s.factor(total);assert total==0
        rows.append(dict(kinds=kind,bridges=0,colors='',split=0,pairings=count,polynomial=str(total)))
    for kinds,k in [('XX',2),('XR',2),('RX',2),('RR',2),('RR',3)]:
        n,m=[2 if c=='X' else 3 for c in kinds];ext=n+m-2*k
        for aa in combinations(range(ext),n-k):
            colors=''.join('A' if i in aa else 'B' for i in range(ext))
            for split in range(ext+1):
                A=list(range(n));BB=list(range(n+split,n+split+m))
                ex=list(range(n,n+split))+list(range(n+split+m,n+m+ext));req=dict(zip(ex,colors))
                total=0;count=0
                for pa in pairings(tuple(range(n+m+ext))):
                    pi={a:b for a,b in pa}|{b:a for a,b in pa}
                    if any(pi[e] not in (A if c=='A' else BB) for e,c in req.items()):continue
                    if any(pi[a] in A for a in A) or any(pi[b] in BB for b in BB):continue
                    if sum(pi[a] in BB for a in A)!=k:continue
                    count+=1;total+=value(pa,kinds,[A,BB])
                total=s.factor(total)
                expected=0
                if kinds=='XX':expected=4*l*(1+q)**2+16*l*l*q*(1+q)
                if kinds=='RR' and k==3:
                    expected=4*l*(1-q*q)*(1+q*q)*(1+q+q*q)**2+16*l*l*q*q*(1+q)*(1+q*q)*(1+q+q*q)*(2*q*q-q+2)
                assert s.expand(total-expected)==0,(kinds,k,colors,split,total)
                rows.append(dict(kinds=kinds,bridges=k,colors=colors,split=split,pairings=count,polynomial=str(total)))
    names=['code/mixed_cubic_clusters.py','code/singlet_coefficient.py','code/singlet_spectral_graphs.py']
    target.write_text(json.dumps(dict(status='passed',templates=len(rows),pairings=sum(r['pairings'] for r in rows),records=rows,
        scope='All singleton and >=2-bridge two-vertex clusters for degree2/3; proof requires external-mark invariance in manuscript',
        environment=dict(python=platform.python_version(),sympy=s.__version__),
        input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names}),indent=2)+'\n')
    print(json.dumps(dict(templates=len(rows),pairings=sum(r['pairings'] for r in rows))))
if __name__=='__main__':main()
