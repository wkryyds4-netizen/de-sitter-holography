"""Pure covariant trace and one-selected-state formulas; no file writes.

H=h, W=sqrt(lambda)u, J=j. Physical gamma=exp(-lambda).
Exact algebraic checks may use independent rational gamma and lambda.
"""
from fractions import Fraction as F
from functools import lru_cache
from bilocal_multisource_chords import compressed_moment,evaluate

@lru_cache(None)
def trace(word,gamma,lam):
    gamma,lam=F(gamma),F(lam)
    return evaluate(compressed_moment(word),gamma**4,lam,word.count('J'))

def selected(word,gamma,lam):
    gamma,lam=F(gamma),F(lam);s=gamma*gamma
    return ((1-s*(1+lam))*trace(word,gamma,lam)
        +s*sum(sign*trace(prefix+word,gamma,lam)
               for prefix,sign in [('HH',1),('HW',-1),('WH',-1),('WW',1)]))

def reduced_hj(word,gamma,lam):
    assert set(word)<=set('HJ')
    gamma,lam=F(gamma),F(lam);s=gamma*gamma
    return (1-s)*trace(word,gamma,lam)+s*trace('HH'+word,gamma,lam)

def curvature_commutator(gamma):
    gamma=F(gamma);q=gamma**4;s=gamma**2;r=q*q
    return 2*(1-r)+s*(1+q)*(1-r)**2

def density_nonstationarity(gamma,lam):
    gamma,lam=F(gamma),F(lam);q=gamma**4
    return 2*gamma**4*lam*(1+lam)*(1+q)*(1-q*q)
