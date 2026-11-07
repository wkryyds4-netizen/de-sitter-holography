"""Exploratory exact diagonal rank-two input at half filling; no old output overwritten."""
import numpy as np
from memory_quadratic_contraction import Contraction
from memory_rank_blocks import casimir_matrix, projectors

o=Contraction(6,2,3)
def diagonal(m):
    return np.diag([(int(bool(x&1))-int(bool(x&2)))*(int(bool(x&4))-int(bool(x&8)))
                    for x in range(64) if x.bit_count()==m]).astype(complex)
h,k=diagonal(3),diagonal(2)
assert np.array_equal(casimir_matrix(6,3)@h.ravel(),20*h.ravel())
pr=projectors(casimir_matrix(6,3),6,3)
for name,x in [('h2',h@h),('A',o.A(h,k)),('Q',o.contract(h,k))]:
    print(name,flush=True)
    for r,(num,den) in enumerate(pr):
        y=num@x.ravel()
        print(r,den,np.max(abs(y)),np.vdot(y,y).real,flush=True)
