"""Finite endpoint-orbit certificate for arbitrary-length cubic two-arc words."""
from pathlib import Path
from itertools import combinations, permutations
from math import factorial
import json,hashlib,platform
import sympy as s
from singlet_spectral_graphs import graph

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-cubic-kernel'
q=s.Symbol('q')
D=(1+q*q)*(1+q+q*q)
B=2*(1-q)*(1+q**3)
weights=[-q*(1+q+q*q),-q*(1+q+q*q),D]
mpairs=[(0,1),(1,2),(0,2)]
def crosses(a,b,c,d):
    a,b=sorted((a,b));c,d=sorted((c,d))
    return int(a<c<b<d or c<a<d<b)
def canon(pa):return tuple(sorted(tuple(sorted(edge)) for edge in pa))
def orbit(k,colors,split):
    """Both contiguous blocks of size 3; colored external endpoints on two arcs."""
    A=list(range(3));Bslots=list(range(3+split,6+split))
    external=list(range(3,3+split))+list(range(6+split,6+len(colors)))
    ae=[e for e,c in zip(external,colors) if c=='A']
    be=[e for e,c in zip(external,colors) if c=='B']
    seed=list(zip(A[:k],Bslots[:k]))+list(zip(A[k:],ae))+list(zip(Bslots[k:],be))
    words=set()
    for ap in permutations(A):
        for bp in permutations(Bslots):
            pm=dict(zip(A,ap))|dict(zip(Bslots,bp))
            words.add(canon((pm.get(a,a),pm.get(b,b)) for a,b in seed))
    assert len(words)==36//factorial(k)
    return A,Bslots,sorted(words)
def factors(pa,A,Bslots):
    pi={a:b for a,b in pa}|{b:a for a,b in pa}
    cr=sum(crosses(a,b,c,d) for i,(a,b) in enumerate(pa) for c,d in pa[i+1:])
    alphas=[]
    for block in [A,Bslots]:
        alphas.append([2-4*crosses(block[a],pi[block[a]],block[b],pi[block[b]]) for a,b in mpairs])
    delta=[[4*int(pi[A[a]]==Bslots[d] and pi[A[b]]==Bslots[c])
      for c,d in mpairs] for a,b in mpairs]
    return cr,alphas,delta
def main():
    OUT.mkdir(exist_ok=True)
    target=OUT/'orbit-certificate.json'
    if target.exists():raise FileExistsError(target)
    records=[];checked=0;words_total=0;failures=[]
    for k in range(4):
        m=6-2*k
        for aset in combinations(range(m),3-k):
            colors=''.join('A' if i in aset else 'B' for i in range(m))
            for split in range(m+1):
                A,Bslots,words=orbit(k,colors,split);words_total+=len(words)
                # Aggregate integer coefficient tensors before symbolic reduction.
                tensor={}
                for wi,pa in enumerate(words):
                    cr,alpha,delta=factors(pa,A,Bslots)
                    slot=tensor.setdefault(cr,[[0]*4 for _ in range(4)])
                    aa=alpha[0]+[1];bb=alpha[1]+[1]
                    for i in range(4):
                        for j in range(4):slot[i][j]+=aa[i]*bb[j]
                    slot_delta=tensor.setdefault(('delta',cr),[[0]*3 for _ in range(3)])
                    for i in range(3):
                        for j in range(3):slot_delta[i][j]+=delta[i][j]
                    # Independent existing determinant-minor rule, representative orbit members.
                    if wi in [0,len(words)-1]:
                        for i,(a,b) in enumerate(mpairs):
                            for j,(c,d) in enumerate(mpairs):
                                g=graph(6+m,pa,((A[a],A[b]),(Bslots[c],Bslots[d])))
                                assert (g['crossings'],g['alpha'],g['beta'],g['delta'])==(cr,alpha[0][i],alpha[1][j],delta[i][j])
                                checked+=1
                ww=weights+[-B]
                quadratic=sum(q**cr*sum(mat[i][j]*ww[i]*ww[j] for i in range(4) for j in range(4))
                  for cr,mat in tensor.items() if isinstance(cr,int))
                linear=sum(q**key[1]*sum(mat[i][j]*weights[i]*weights[j] for i in range(3) for j in range(3))
                  for key,mat in tensor.items() if isinstance(key,tuple))
                quadratic=s.factor(quadratic);linear=s.factor(linear)
                if k<3 and (quadratic!=0 or linear!=0):failures.append((k,colors,split,str(linear),str(quadratic)))
                records.append(dict(bridges=k,external_colors=colors,arc_split=split,orbit_size=len(words),
                  lambda_coefficient_numerator=str(linear),lambda_squared_coefficient_numerator=str(quadratic)))
        print('bridges',k,'done; failures',len(failures),flush=True)
    data=dict(status='passed' if not failures else 'counterexamples',denominator=str(s.expand(D**2)),
      templates=len(records),orbit_words=words_total,local_minor_comparisons=checked,records=records,failures=failures,
      environment=dict(python=platform.python_version(),sympy=s.__version__),
      input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in ['code/cubic_kernel_orbits.py','code/singlet_spectral_graphs.py','code/singlet_coefficient.py']})
    target.write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps({key:data[key] for key in ['status','templates','orbit_words','local_minor_comparisons','failures']}),flush=True)
    print('surviving',records[-1],flush=True)
if __name__=='__main__':main()
