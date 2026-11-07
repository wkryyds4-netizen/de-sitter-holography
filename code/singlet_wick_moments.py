"""Finite-N Gaussian contractions for the same-sample derivative singlet.

All matrix algebra uses Gaussian integers; integer bounds are checked before
norm sums. Fraction arithmetic restores covariance and trace denominators.
No disorder sampling and no spectral or charge cutoff is used.
"""
from fractions import Fraction
from math import comb
import numpy as np
from complex_syk_dictionary import fermions, creators, covariance, a2_operator


def gue_basis(n, p):
    cs = fermions(n)
    bs = creators(cs, p)
    vs, r = [], []
    for a in range(len(bs)):
        vs.append(bs[a] @ bs[a].conj().T)
        r.append(2)
        for b in range(a):
            t = bs[a] @ bs[b].conj().T
            vs.extend([t + t.conj().T, 1j * (t - t.conj().T)])
            r.extend([1, 1])
    return cs, np.array(vs), np.array(r, dtype=np.int64)


def integer_norm_by_charge(a, charges, n):
    """Return sum of |entries|^2 in each input-charge column (unnormalized)."""
    assert np.all(a.real == np.rint(a.real))
    assert np.all(a.imag == np.rint(a.imag))
    # Conservative bound: also keeps float matmul results far below 2**53.
    assert np.max(np.abs(a)) < 2**20
    rr, ii = a.real.astype(np.int64), a.imag.astype(np.int64)
    values = rr * rr + ii * ii
    assert values.size * 2**41 < 2**63
    return np.array([values[..., charges == m].sum(dtype=np.int64)
                     for m in range(n + 1)], dtype=np.int64)


def compute(n, p):
    cs, v, r = gue_basis(n, p)
    count, d, _ = v.shape
    charges = np.array([s.bit_count() for s in range(d)])
    deriv = np.array([[h @ c - c @ h for c in cs] for h in v])
    raw = np.einsum('aiyx,biyz->abxz', deriv.conj(), deriv)
    b2 = raw + raw.swapaxes(0, 1)  # 2 B_ab, all Gaussian integers
    mean4 = np.einsum('a,aaxy->xy', r, b2)  # 4 E[A]
    bsum = np.zeros(n + 1, dtype=np.int64)
    for a in range(count):
        for b in range(count):
            bsum += int(r[a] * r[b]) * integer_norm_by_charge(b2[a,b], charges, n)
    csum = np.zeros(n + 1, dtype=np.int64)
    contractions = np.zeros_like(v)  # sum_a r_a (6 C_aac)
    # Ordered triple sum, batched only over the last index.
    for a in range(count):
        for b in range(count):
            c6 = (v[a] @ b2[b] - b2[b] @ v[a]
                  + v[b] @ b2[a] - b2[a] @ v[b]
                  + v @ b2[a,b] - b2[a,b] @ v)
            for c in range(count):
                csum += int(r[a] * r[b] * r[c]) * integer_norm_by_charge(c6[c], charges, n)
            if a == b:
                contractions += int(r[a]) * c6
    lsum = np.zeros(n + 1, dtype=np.int64)
    for c in range(count):
        lsum += int(r[c]) * integer_norm_by_charge(contractions[c], charges, n)
    meansq = integer_norm_by_charge(mean4, charges, n)
    rows = []
    for m in range(n + 1):
        dm = comb(n,m)
        mean = Fraction(int(np.trace(mean4[np.ix_(charges==m,charges==m)]).real),4*dm)
        mu0 = Fraction(int(meansq[m]),16*dm) + Fraction(int(bsum[m]),8*dm)
        mu2c = Fraction(int(csum[m]),48*dm)
        mu2l = Fraction(int(lsum[m]),32*dm)
        rows.append(dict(m=m, dimension=dm, mean=mean, mu0=mu0,
                         ensemble_centered_gram=mu0-mean**2,
                         mu2_cubic=mu2c, mu2_linear=mu2l, mu2=mu2c+mu2l))
    total = {key:sum(Fraction(row['dimension'],d)*row[key] for row in rows)
             for key in ['mean','mu0','ensemble_centered_gram','mu2_cubic','mu2_linear','mu2']}
    # Deterministic polynomial check against original commutator definition.
    x = np.array([((7*a+3)%13-6)/7 for a in range(count)])
    h = np.einsum('a,axy->xy',x,v)
    aa = a2_operator(h,cs)
    poly = np.einsum('a,b,abxy->xy',x,x,b2)/2
    residual = np.linalg.norm(aa-poly)/max(1,np.linalg.norm(aa))
    return dict(N=n,p=p,real_coordinates=count,sector=rows,total=total,
                polynomial_relative_residual=float(residual),
                integer_contraction_max=float(np.max(np.abs(contractions)))), (cs,v,r,b2)


def independent_wick6(v, r, b2):
    """All 15 pairings of unsymmetrized D_abc=[V_a,B_bc], small controls.

    Independent of the 6+9 symmetric-tensor decomposition in compute().
    """
    def pairs(indices):
        if not indices:
            yield []
        else:
            a=indices[0]
            for j in range(1,len(indices)):
                for rest in pairs(indices[1:j]+indices[j+1:]):
                    yield [(a,indices[j])]+rest
    count,d,_=v.shape
    d2 = np.array([h@b2-b2@h for h in v])
    terms=[]
    for pairing in pairs(tuple(range(6))):
        val=0.0
        for a in range(count):
            for b in range(count):
                for c in range(count):
                    slots=[None]*6
                    for (i,j),idx in zip(pairing,(a,b,c)):
                        slots[i]=slots[j]=idx
                    lhs=d2[tuple(slots[:3])]
                    rhs=d2[tuple(slots[3:])]
                    val += int(r[a]*r[b]*r[c])*np.vdot(lhs,rhs).real
        assert val == round(val) and abs(val)<2**53
        terms.append(Fraction(int(val),32*d))
    return terms
