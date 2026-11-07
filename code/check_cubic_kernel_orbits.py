"""Independent enumeration: all perfect matchings and integer polynomials."""
from pathlib import Path
from math import factorial
import json,hashlib
from singlet_coefficient import pairings
from singlet_spectral_graphs import graph
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-cubic-kernel'
# Ascending integer polynomial coefficients; no symbolic factoring in this verifier.
W=[(0,-1,-1,-1),(0,-1,-1,-1),(1,1,2,1,1)]
C=(-2,2,0,-2,2)
def mul(a,b):
    r=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):r[i+j]+=x*y
    return r
def add(dest,source,scale=1,shift=0):
    for i,x in enumerate(source):dest[i+shift]+=scale*x
def main():
    target=OUT/'independent-orbits.json'
    if target.exists():raise FileExistsError(target)
    groups={};total=0;minor_calls=0
    for k in range(4):
        m=6-2*k
        for split in range(m+1):
            length=6+m;A=set(range(3));B=set(range(3+split,6+split))
            ext=sorted(set(range(length))-A-B)
            am=[(0,1),(1,2),(0,2)]
            bm=[tuple(x+3+split for x in pair) for pair in am]
            for pa in pairings(tuple(range(length))):
                pi={a:b for a,b in pa}|{b:a for a,b in pa}
                if any(pi[a] in A for a in A) or any(pi[b] in B for b in B):continue
                if any(pi[e] in ext for e in ext):continue
                if sum(pi[a] in B for a in A)!=k:continue
                colors=''.join('A' if pi[e] in A else 'B' for e in ext)
                key=(k,colors,split);row=groups.setdefault(key,dict(words=0,linear=[0]*30,quadratic=[0]*30))
                row['words']+=1;total+=1
                fa=list(C);fb=list(C);delta=[0]*9
                for i,ma in enumerate(am):
                    for j,mb in enumerate(bm):
                        g=graph(length,pa,(ma,mb));minor_calls+=1
                        if j==0:add(fa,W[i],g['alpha'])
                        if i==0:add(fb,W[j],g['beta'])
                        add(delta,mul(W[i],W[j]),g['delta'])
                add(row['linear'],delta,shift=g['crossings'])
                add(row['quadratic'],mul(fa,fb),shift=g['crossings'])
        print('independent matchings k',k,'complete',flush=True)
    old=json.loads((OUT/'orbit-certificate.json').read_text())
    assert len(groups)==old['templates']==177 and total==old['orbit_words']==6234
    for key,row in groups.items():
        assert row['words']==36//factorial(key[0])
        if key[0]<3:assert not any(row['linear']+row['quadratic']),key
    last=groups[(3,'',0)]
    # Verify surviving coefficients directly against expanded integer products.
    expected_l=mul(mul((4,0,-4),(1,0,1)),mul((1,1,1),(1,1,1)))
    expected_q=mul(mul(mul((0,0,16,16),(1,0,1)),(1,1,1)),(2,-1,2))
    assert last['linear'][:len(expected_l)]==expected_l and not any(last['linear'][len(expected_l):])
    assert last['quadratic'][:len(expected_q)]==expected_q and not any(last['quadratic'][len(expected_q):])
    result=dict(status='passed',method='all perfect matchings, determinant-minor rule and integer coefficient arrays; no orbit generator or symbolic factorizer imported',
      templates=len(groups),matching_words=total,minor_calls=minor_calls,cancelled_templates=176,
      nonzero_linear_coefficients=expected_l,nonzero_quadratic_coefficients=expected_q,
      input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in ['code/check_cubic_kernel_orbits.py','code/singlet_coefficient.py','code/singlet_spectral_graphs.py','results/susskind-cubic-kernel/orbit-certificate.json']})
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({key:result[key] for key in ['status','templates','matching_words','minor_calls','cancelled_templates']}))
if __name__=='__main__':main()
