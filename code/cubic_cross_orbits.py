"""All-order R-S and R-Z cancellations on the existing exhaustive templates."""
from pathlib import Path
from itertools import combinations
from collections import defaultdict
import hashlib,json
import sympy as s
from cubic_kernel_orbits import orbit,factors,weights,B,q,mpairs
from singlet_spectral_graphs import graph
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-cubic-cross'
def main():
    target=OUT/'cross-orbits.json'
    if target.exists():raise FileExistsError(target)
    rows=[];total=0;minor_checks=0
    for k in range(4):
        m=6-2*k
        for aset in combinations(range(m),3-k):
            colors=''.join('A' if i in aset else 'B' for i in range(m))
            for split in range(m+1):
                A,BB,words=orbit(k,colors,split);total+=len(words)
                ll=defaultdict(lambda:[0]*3);qq=defaultdict(lambda:[0]*4);zz=defaultdict(lambda:[0]*4)
                for pa in words:
                    cr,al,de=factors(pa,A,BB)
                    for i,(a,b) in enumerate(mpairs):
                        for j,(c,d) in enumerate(mpairs[:2]):
                            g=graph(6+m,pa,((A[a],A[b]),(BB[c],BB[d])))
                            assert (g['crossings'],g['alpha'],g['beta'],g['delta'])==(cr,al[0][i],al[1][j],de[i][j]);minor_checks+=1
                    for i in range(3):ll[cr][i]+=de[i][0]+de[i][1]
                    for i,aa in enumerate(al[0]+[1]):
                        qq[cr][i]+=aa*(al[1][0]+al[1][1]);zz[cr][i]+=aa
                ww=weights+[-B]
                lin=s.expand(sum(q**cr*sum(v[i]*weights[i] for i in range(3)) for cr,v in ll.items()))
                quad=s.expand(sum(q**cr*sum(v[i]*ww[i] for i in range(4)) for cr,v in qq.items()))
                zzv=s.expand(sum(q**cr*sum(v[i]*ww[i] for i in range(4)) for cr,v in zz.items()))
                assert lin==quad==zzv==0
                rows.append(dict(bridges=k,colors=colors,split=split,orbit_size=len(words),RS_lambda=str(lin),RS_lambda2=str(quad),RZ_lambda=str(zzv)))
        print('cross orbit bridges',k,'passed',flush=True)
    names=['code/cubic_cross_orbits.py','code/cubic_kernel_orbits.py','code/singlet_spectral_graphs.py','code/singlet_coefficient.py']
    target.write_text(json.dumps(dict(status='passed',templates=len(rows),words=total,minor_checks=minor_checks,records=rows,
      scope='Same exhaustive spectator-independent templates as cubic-kernel; changed bilinear polynomial weights to R-S and R-Z. All three coefficient arrays vanish.',
      input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names}),indent=2)+'\n')
if __name__=='__main__':main()
