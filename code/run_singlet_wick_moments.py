"""One-shot exact covariance sums and independent degree-six checks."""
from pathlib import Path
from fractions import Fraction
from itertools import product
from math import factorial, pi
import hashlib
import json
import platform
import numpy as np
from singlet_wick_moments import compute, independent_wick6, gue_basis
from complex_syk_dictionary import a2_operator, charge_spectrum

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-singlet-moments'


def check_subspace():
    cs, full, r=gue_basis(4,2)
    indices=[2,10,21]
    v=full[indices]*np.sqrt(r[indices]/2)[:,None,None]
    deriv=np.array([[h@c-c@h for c in cs] for h in v])
    raw=np.einsum('aiyx,biyz->abxz',deriv.conj(),deriv)
    b=(raw+raw.swapaxes(0,1))/2
    c=np.array([[(v[a]@b[z]-b[z]@v[a]+v[z]@b[a]-b[a]@v[z]
                 +v@b[a,z]-b[a,z]@v)/3 for z in range(3)] for a in range(3)])
    mean=np.einsum('aaxy->xy',b)
    first=np.einsum('aabxy->bxy',c)
    mu0=(np.vdot(mean,mean).real+2*np.vdot(b,b).real)/16
    cubic=6*np.vdot(c,c).real/16
    linear=9*np.vdot(first,first).real/16
    nodes,weights=np.polynomial.hermite_e.hermegauss(4)
    weights=weights/np.sqrt(2*pi)
    q0=q2=0.0
    for inds in product(range(4),repeat=3):
        x=nodes[list(inds)];weight=float(np.prod(weights[list(inds)]))
        h=np.einsum('a,axy->xy',x,v)
        aa=a2_operator(h,cs);kk=h@aa-aa@h
        q0+=weight*np.vdot(aa,aa).real/16
        q2+=weight*np.vdot(kk,kk).real/16
    # One deterministic realization also checks the frequency-moment relation.
    h=np.einsum('a,axy->xy',[.3,-.8,1.2],v)
    aa=a2_operator(h,cs);kk=h@aa-aa@h
    omega,weight,_,_=charge_spectrum(h,aa,4)
    direct=float(np.vdot(kk,kk).real/16)
    spectral=float(np.dot(omega**2,weight))
    residuals=[abs(q0-mu0)/max(1,abs(mu0)),abs(q2-cubic-linear)/max(1,abs(cubic+linear)),
               abs(direct-spectral)/max(1,abs(direct))]
    assert max(residuals)<1e-12
    return dict(indices=indices,gaussian_nodes=64,mu0_tensor=float(mu0),mu0_quadrature=q0,
                mu2_cubic=float(cubic),mu2_linear=float(linear),mu2_quadrature=q2,
                direct_frequency_moment=direct,spectral_frequency_moment=spectral,
                scaled_residuals=residuals)


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    target=OUT/'moments.json'
    if target.exists(): raise FileExistsError(target)
    records=[]
    for n,p in [(3,1),(3,2),(4,2)]:
        record,components=compute(n,p)
        w=4**p*Fraction(n,p*p)*Fraction(factorial(p),n**p)**2
        record['standard_CAR_coefficient_variance_Js1']=w
        record['string_mu0_Js1']=record['total']['mu0']*w**2
        record['string_mu2_Js1']=record['total']['mu2']*w**3
        assert record['polynomial_relative_residual']<1e-12
        if p==1 or n<=p+1: assert record['total']['mu2']==0
        if (n,p)==(4,2):
            terms=independent_wick6(components[1],components[2],components[3])
            assert sum(terms)==record['total']['mu2']
            record['independent_15_pairings']=terms
            record['pairing_cancellation_ratio']=sum(abs(t) for t in terms)/sum(terms)
        records.append(record)
    payload=dict(status='passed',method='Complete finite Gaussian Wick sums; no random sampling',
                 conventions='Raw tables have w=1; trace 2^-N Tr. Sector tables use binomial-normalized traces. Total ensemble_centered_gram subtracts E[A]|Q, not sample trace or full conserved projection.',
                 environment=dict(python=platform.python_version(),numpy=np.__version__),
                 records=records,independent_cubature=check_subspace(),
                 input_sha256={str(p.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest()
                     for p in [ROOT/'code/singlet_wick_moments.py',Path(__file__),ROOT/'code/complex_syk_dictionary.py']})
    with target.open('x',encoding='utf-8') as f: json.dump(payload,f,indent=2,default=str);f.write('\n')
    print(json.dumps(dict(status='passed',mu2=str(records[-1]['string_mu2_Js1']),
                         cubature=payload['independent_cubature']),indent=2))


if __name__=='__main__':main()
