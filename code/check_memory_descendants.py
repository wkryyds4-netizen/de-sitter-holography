"""Independent finite CAR checks of raw descendants and their single contractions."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,platform,time
import numpy as np
import sympy as s
from memory_quadratic_contraction import Contraction,adj,anti
from check_memory_two_charge_filter import words

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-descendants'


def raw(o,h,k):
    # Direct body-rank and derivative definitions, not a five-word coefficient table.
    h2=h@h;h3=h2@h;h4=h2@h2
    a=o.A(h,k);t3=o.rem(k@k@k)+(3*o.p-o.m)*h3
    t4=o.rem(k@k@k@k)+(4*o.p-o.m)*h4
    E=anti(a,h2)/2;M=h@a@h;S=anti(h,t3)/2-E-M
    C=t4-2*E-M-2*S
    return np.array([E,M,S,C,h4])


def contraction(o,h,k):
    d=words(o,h,k);p,m=o.p,o.m
    return np.array([d['D22']/2+(2*p-m)*d['D40'],
                     d['D121']+(2*p-m)*d['D40'],
                     d['D13']/2-d['D22']/2-d['D121']+(m-p)*d['D40'],
                     d['D4']-d['D13']+d['D121'],d['D40']])


def integer_check(z):
    assert np.max(abs(z))<2**45
    assert np.array_equal(z.real,np.rint(z.real))
    assert np.array_equal(z.imag,np.rint(z.imag))


def independent(o,h,k,base_single):
    base=raw(o,h,k);value=np.zeros_like(base)
    for rv,hh,kk,single in zip(o.rs,o.hs,o.ks,base_single):
        delta=raw(o,h+hh,k+kk)+raw(o,h-hh,k-kk)-2*base-2*single
        integer_check(2*delta)
        value+=rv*delta/4
    return value


def main():
    target=OUT/'finite-descendants.json'
    if target.exists():raise FileExistsError(target)
    rows=[];start=time.time()
    for N,p,m,full in [(3,1,2,True),(4,2,3,True),(6,2,3,False)]:
        o=Contraction(N,p,m);nr=len(o.rs)
        singles=[raw(o,h,k) for h,k in zip(o.hs,o.ks)]
        directions=[(o.hs[a],o.ks[a]) for a in range(nr)] if full else []
        if full:
            directions += [(o.hs[a]+o.hs[b],o.ks[a]+o.ks[b]) for a in range(nr) for b in range(a+1,nr)]
        for seed in [3,7,11]:
            x=np.array([((seed+2)*i+seed)%9-4 for i in range(nr)])
            directions.append((np.einsum('a,aij->ij',x,o.hs),np.einsum('a,aij->ij',x,o.ks)))
        maxentry=0
        for h,k in directions:
            value=contraction(o,h,k);ind=independent(o,h,k,singles)
            integer_check(8*value);integer_check(8*ind)
            assert np.array_equal(value,ind),(N,p,m)
            rr=raw(o,h,k)
            assert np.array_equal(rr[3]-2*rr[2]+(rr[0]+rr[1])/2,o.core(h,k))
            maxentry=max(maxentry,float(np.max(abs(8*value))))
        # A second contraction recovers the scalar mean exactly in every coordinate.
        mean=sum(rv*contraction(o,h,k)/4 for rv,h,k in zip(o.rs,o.hs,o.ks))
        scalar=[F(int(np.trace(z).real*8),8*o.d) for z in mean]
        for z,sc in zip(mean,scalar):assert np.array_equal(z,float(sc)*np.eye(o.d))
        row=dict(N=N,p=p,m=m,coordinate_count=nr,full_quadratic_tensor=full,
                 directions=len(directions),contraction_identities=5*len(directions),
                 core_raw_identities=len(directions),means=list(map(str,scalar)),max_scaled_entry=maxentry)
        rows.append(row);print(row,flush=True)
    names=['code/check_memory_descendants.py','code/memory_quadratic_contraction.py',
           'code/check_memory_two_charge_filter.py','code/singlet_wick_moments.py']
    data=dict(status='passed',rows=rows,elapsed_seconds=time.time()-start,
              method='Quartic even polarization over every independent Gaussian direction, using direct descendant definitions. Complete quadratic coefficient tensor in the first two systems; three deterministic larger half-filled directions. Checked dyadic entries remain within exact integer range.',
              environment=dict(python=platform.python_version(),numpy=np.__version__,sympy=s.__version__),
              input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names})
    target.write_text(json.dumps(data,indent=2)+'\n')


if __name__=='__main__':main()
