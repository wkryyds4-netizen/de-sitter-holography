"""Independent fixed-bubble chord identities and conserved return quadrature.

Finite polynomial checks do not substitute for the same-J tag argument.
No historical drivers are run; fresh output only.
"""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import factorial
import json,platform
import numpy as np
import scipy
from scipy.linalg import eigh_tridiagonal
from numpy.polynomial.legendre import leggauss

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-return-background'

@lru_cache(None)
def pairings(points):
    if not points: return ((),)
    first=points[0];out=[]
    for k in range(1,len(points)):
        for rest in pairings(points[1:k]+points[k+1:]):
            out.append(((first,points[k]),)+rest)
    return tuple(out)

def crossing(a,b):
    x,y=a;u,v=b
    return int(x<u<y<v or u<x<v<y)

def p_add(a,b,scale=F(1)):
    out=list(a)+[F(0)]*max(0,len(b)-len(a))
    for k,x in enumerate(b):out[k]+=scale*x
    return out

def p_mul(a,b):
    out=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    return out

def conditional_monomials(q,maximum):
    h=[[F(1)],[F(0),F(1)]]
    for n in range(1,maximum):
        number=sum(q**j for j in range(n))
        h.append(p_add([F(0)]+h[-1],h[-2],-number))
    out=[]
    for n in range(maximum+1):
        rem=[F(0)]*n+[F(1)]; image=[F(0)]*(n+1)
        for k in range(n,-1,-1):
            coeff=rem[k]
            image=p_add(image,h[k],coeff*q**k)
            rem=p_add(rem,h[k],-coeff)
        assert all(x==0 for x in rem)
        out.append(image)
    return out

@lru_cache(None)
def fixed_bubble_counts(exponents):
    a,b,c,d,e=exponents; pos=0
    spectators=list(range(a));pos=a;one=pos;pos+=1
    spectators+=list(range(pos,pos+b));pos+=b;two=pos;pos+=1
    spectators+=list(range(pos,pos+c));pos+=c;three=pos;pos+=1
    spectators+=list(range(pos,pos+d));pos+=d;four=pos;pos+=1
    spectators+=list(range(pos,pos+e));chords=((one,two),(three,four))
    counts={}
    for pp in pairings(tuple(spectators)):
        k=sum(crossing(x,y) for j,x in enumerate(pp) for y in pp[j+1:])
        k+=sum(crossing(x,y) for x in pp for y in chords)
        counts[k]=counts.get(k,0)+1
    return counts

def chord_checks():
    rows=[]; degrees=8
    configs=[x for x in product(range(degrees+1),repeat=5) if sum(x)<=degrees and sum(x)%2==0]
    for q in [F(1,4),F(1,2),F(3,4)]:
        image=conditional_monomials(q,degrees)
        moments=[]
        for n in range(degrees+1):
            val=F(0)
            if n%2==0:
                for pp in pairings(tuple(range(n))):
                    k=sum(crossing(x,y) for j,x in enumerate(pp) for y in pp[j+1:])
                    val+=q**k
            moments.append(val)
        for exponents in configs:
            a,b,c,d,e=exponents
            counts=fixed_bubble_counts(exponents)
            chords=sum((count*q**k for k,count in counts.items()),F(0))
            poly=[F(0)]*(a+c+e)+p_mul(image[b],image[d])
            expected=sum((x*moments[j] for j,x in enumerate(poly)),F(0))
            assert chords==expected,(q,exponents,chords,expected)
        rows.append(dict(q=str(q),source_word_configurations=len(configs),max_H_degree=degrees))
    # This cross contraction would survive for one common Gaussian probe.
    # It is precisely the term the microscopic tag-intersection proof must suppress.
    return dict(rows=rows,total_exact_equalities=len(configs)*3,
                all_four_slots_same_probe_moment='2+q; fixed two internal bubbles contribute only 1',
                note='Tests q-conditional identity for prescribed chords, not the same-J DS counting itself.')

def q_basis(x,q,sigma,count):
    out=np.zeros((count,len(np.atleast_1d(x))));out[0]=1
    if count==1:return out
    beta=lambda n:sigma*np.sqrt((1-q**n)/(1-q))
    out[1]=x/beta(1)
    for n in range(1,count-1):out[n+1]=(x*out[n]-beta(n)*out[n-1])/beta(n+1)
    return out

def scalar_return(lam,energy_fraction,s,order,quad):
    q=np.exp(-4*lam);sigma=1/np.sqrt(lam);edge=2*sigma/np.sqrt(1-q)
    e,z=eigh_tridiagonal(np.zeros(order),sigma*np.sqrt((1-q**np.arange(1,order))/(1-q)))
    weights=z[0]**2;target=np.array([energy_fraction*edge])
    basis_nodes=q_basis(e,q,sigma,order);basis_target=q_basis(target,q,sigma,order)[:,0]
    # Mehler conditional density at the target energy, integrated against nu_q.
    conditional=weights*((basis_target*q**np.arange(order))@basis_nodes)
    assert abs(conditional.sum()-1)<5e-13
    assert conditional.min()>-5e-13
    def g(times):
        times=np.asarray(times)
        return np.exp(1j*times[:,None]*(e[None,:]-target[0]))@conditional
    nodes,w=leggauss(quad);v=(nodes+1)*abs(s)/2;wv=w*abs(s)/2
    integ=0j
    for vi,wi in zip(v,wv):
        ww=(nodes+1)*(abs(s)-vi)/2;weights_w=w*(abs(s)-vi)/2
        integ+=wi*g([vi])[0]*np.dot(weights_w,g(ww))
    value=8/lam*integ
    if s<0:value=value.conjugate()
    return value

def scalar_checks():
    rows=[]
    for lam,fraction,s in product([.25,.5,1.],[0.,.4,.85],[.5,1.,2.]):
        coarse=scalar_return(lam,fraction,s,64,20)
        fine=scalar_return(lam,fraction,s,96,32)
        err=abs(coarse-fine);assert err<3e-11
        rows.append(dict(lambda_value=lam,Js=1.,energy_fraction_of_edge=fraction,s=s,
                         kappa=[float(fine.real),float(fine.imag)],absolute_order_difference=float(err),
                         conservative_unitarity_bound=4*s*s/lam))
        assert abs(fine)<=4*s*s/lam+1e-12
    return rows

if __name__=='__main__':
    OUT.mkdir(parents=True,exist_ok=True);dest=OUT/'bubble-and-return-checks.json'
    if dest.exists():raise FileExistsError(dest)
    data=dict(chords=chord_checks(),return_background=scalar_checks(),
              environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
              parameters=dict(q='exp(-4 lambda)',variance='Js^2/lambda',spectral_orders=[64,96],simplex_orders=[20,32]),
              scope='Exact rational chord identities and float64 conserved-background quadrature; no finite-N asymptotic fit or certified quadrature error.')
    dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status='passed',exact_equalities=data['chords']['total_exact_equalities'],
                         scalar_values=len(data['return_background']),max_order_difference=max(x['absolute_order_difference'] for x in data['return_background']))))
