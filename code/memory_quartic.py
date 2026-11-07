"""Second centered memory coefficient and quartic single-mark Gram."""
from pathlib import Path
from itertools import combinations
from fractions import Fraction
import json,hashlib,platform
import sympy as s
import numpy as np
from complex_syk_dictionary import sample_hamiltonian,comm
from bilocal_channel import blocks
from memory_force_bridge import exact
from singlet_coefficient import pairings
from singlet_spectral_graphs import graph

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-quartic'
q,lam=s.symbols('q lambda',positive=True)
pairs=list(combinations(range(4),2))
weights=[{(0,1):Fraction(1,2),(2,3):Fraction(1,2)},
         {(1,2):Fraction(1)},
         {(0,2):Fraction(1,2),(1,3):Fraction(1,2)},
         {(0,3):Fraction(1)},{}]


def exact_cases():
    rows=[]
    for N,p,m in [(3,1,2),(4,2,2),(5,2,3)]:
        cs,J,H=sample_hamiltonian(N,p,110660+N)
        # Quantize Hermitian matrix entries to exact Gaussian integers while
        # retaining a p-body construction, independently of the Gaussian draw.
        from complex_syk_dictionary import creators
        bs=creators(cs,p);J=np.rint(J.real)+1j*np.rint(J.imag)
        H=sum(J[a,b]*(aa@bb.conj().T) for a,aa in enumerate(bs) for b,bb in enumerate(bs))
        inds,h,C,K=blocks(N,p,H,cs,m);h,C,K=map(exact,[h,C,K]);I=s.eye(h.rows)
        f=[[c] for c in cs]
        for j in range(4):
            for row in f:row.append(comm(H,row[-1]))
        BB=lambda a,b:exact(sum(row[a].conj().T@row[b] for row in f)[np.ix_(inds,inds)])
        A=BB(1,1);T3=(BB(0,3)+BB(3,0))/2;T4=(BB(0,4)+BB(4,0))/2;B22=BB(2,2)
        LH=lambda x:h*x-x*h
        JH=lambda x:(h*x+x*h)/2
        r=s.Rational(p,m)
        B={0:I,**{j:C.conjugate().T*K**j*C/m for j in range(1,5)}}
        beta={}
        for n in range(1,5):beta[n]=B[n]-sum((beta[j]*B[n-j] for j in range(1,n)),s.zeros(h.rows))
        G=beta[2];Z2=beta[4]-2*JH(beta[3])+JH(JH(G))
        U4=m*(B[4]-4*JH(B[3])+6*JH(JH(B[2]))-4*JH(JH(JH(B[1])))+h**4)
        formula=B22+(1+r)**2*LH(LH(A))/4+2*r*JH(T3)+3*r*r*JH(JH(A))-A*A/m-m*r**4*h**4
        residuals=[T4-B22-2*LH(LH(A)),BB(0,4)-T4-2*LH(T3),
                   U4-B22-LH(LH(A))/2,m*Z2-formula,
                   Z2-Z2.conjugate().T]
        assert all(all(s.expand(v)==0 for v in x) for x in residuals)
        rows.append(dict(N=N,p=p,m=m,identities=5,mZ2_norm_squared=str(s.trace((m*Z2)**2)/h.rows)))
    return rows


def crossed(a,b,c,d):
    a,b=sorted((a,b));c,d=sorted((c,d))
    return int(a<c<b<d or c<a<d<b)


def gram(gap=0,end=0):
    length=8+gap+end;left=list(range(4));right=list(range(4+gap,8+gap))
    out={};count=minor=0
    for pa in pairings(tuple(range(length))):
        if any((a in left and b in left) or (a in right and b in right) for a,b in pa):continue
        count+=1;pi={a:b for a,b in pa}|{b:a for a,b in pa}
        cr=sum(crossed(a,b,c,d) for j,(a,b) in enumerate(pa) for c,d in pa[j+1:])
        alpha=[]
        for block in [left,right]:
            alpha.append([sum(w*(2-4*crossed(block[a],pi[block[a]],block[b],pi[block[b]])) for (a,b),w in ww.items()) for ww in weights])
        tensor=out.setdefault(cr,[[[Fraction(0) for _ in range(5)] for _ in range(5)] for _ in range(3)])
        for i in range(5):
            for j in range(5):
                if i==4 and j==4:tensor[0][i][j]+=1
                elif i==4:tensor[1][i][j]+=alpha[1][j]
                elif j==4:tensor[1][i][j]+=alpha[0][i]
                else:
                    delta=sum(4*wi*wj*int(pi[left[a]]==right[d] and pi[left[b]]==right[c])
                              for (a,b),wi in weights[i].items() for (c,d),wj in weights[j].items())
                    tensor[1][i][j]+=delta;tensor[2][i][j]+=alpha[0][i]*alpha[1][j]
        if count<=3:
            for ma in pairs:
                for mb in pairs:
                    g=graph(length,pa,(tuple(left[i] for i in ma),tuple(right[i] for i in mb)))
                    aa=2-4*crossed(left[ma[0]],pi[left[ma[0]]],left[ma[1]],pi[left[ma[1]]])
                    bb=2-4*crossed(right[mb[0]],pi[right[mb[0]]],right[mb[1]],pi[right[mb[1]]])
                    dd=4*int(pi[left[ma[0]]]==right[mb[1]] and pi[left[ma[1]]]==right[mb[0]])
                    assert (g['crossings'],g['alpha'],g['beta'],g['delta'])==(cr,aa,bb,dd)
                    minor+=1
    mat=s.Matrix(5,5,lambda i,j:s.factor(sum(q**cr*sum(lam**k*s.Rational(t[k][i][j]) for k in range(3)) for cr,t in out.items())))
    return mat,count,minor


def main():
    OUT.mkdir(parents=True,exist_ok=True);target=OUT/'quartic.json'
    if target.exists():raise FileExistsError(target)
    finite=exact_cases();print('Exact second memory identities passed',flush=True)
    G,count,minor=gram();assert G==G.T
    (OUT/'gram-intermediate.json').write_text(json.dumps(dict(gram=[[str(v) for v in row] for row in G.tolist()],pairings=count),indent=2)+'\n')
    print('Gram built; determinant',s.factor(G.det()),flush=True)
    nuisance=[0,1,2,4]
    A=G.extract(nuisance,nuisance);b=G.extract(nuisance,[3])
    coeff=A.inv()*b;coeff=coeff.applyfunc(s.factor)
    v=s.zeros(5,1);v[3]=1
    for j,i in enumerate(nuisance):v[i]=-coeff[j]
    var=s.factor((v.T*G*v)[0]);print('Projection',list(coeff),'variance',var,flush=True)
    G20,n20,m20=gram(0,2);G11,n11,m11=gram(1,1)
    K=G20+G20.T-2*G11
    kin=s.factor((v.T*K*v)[0]);mu=s.factor(kin/var/lam)
    print('mu2/Js2',mu,flush=True)
    names=['code/memory_quartic.py','code/complex_syk_dictionary.py','code/bilocal_channel.py',
           'code/memory_force_bridge.py','code/singlet_coefficient.py','code/singlet_spectral_graphs.py']
    result=dict(status='passed',finite=finite,basis=['edge=(01+23)/2','middle=12','skip=(02+13)/2','long=03','plain'],
                gram=[[str(vv) for vv in row] for row in G.tolist()],projection=list(map(str,coeff)),
                residual_vector=list(map(str,v)),residual_variance=str(var),residual_kinetic=str(kin),
                normalized_mu2_over_Js2=str(mu),pairings=[count,n20,n11],minor_comparisons=minor+m20+m11,
                environment=dict(python=platform.python_version(),numpy=np.__version__,sympy=s.__version__),
                input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names})
    target.write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
