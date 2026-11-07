"""Finite Gaussian-chaos projection of raw T3; exploratory exact-coordinate data."""
from pathlib import Path
import json,hashlib
import numpy as np
from singlet_wick_moments import gue_basis
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-cubic-linear'
def main():
    OUT.mkdir(exist_ok=True)
    target=OUT/'finite-probe.json'
    if target.exists():raise FileExistsError(target)
    rows=[]
    for n,p in [(3,1),(4,2),(5,2)]:
        cs,vs,rs=gue_basis(n,p);dim=1<<n
        Q=np.diag([a.bit_count() for a in range(dim)])
        def D(Y):return Q@Y-sum(c.conj().T@Y@c for c in cs)
        def twiceT(H):
            A=2*p*H@H-D(H@H)
            return 6*p*H@H@H-2*D(H@H@H)-3*(H@A+A@H)
        x=np.array([(7*i+3)%7-3 for i in range(len(vs))])
        # Nonuniform deterministic integer coefficients.
        x=np.array([(7*i+3)%11-5 for i in range(len(vs))])
        H=np.einsum('a,aij->ij',x,vs);t0=twiceT(H)
        linear8=sum(int(r)*(twiceT(H+v)+twiceT(H-v)-2*t0) for v,r in zip(vs,rs))
        linear=linear8/8
        comm=H@linear-linear@H
        for m in range(n+1):
            ids=[a for a in range(dim) if a.bit_count()==m]
            h=H[np.ix_(ids,ids)];f=linear[np.ix_(ids,ids)]
            design=np.stack([h.ravel(),np.eye(len(ids)).ravel()],axis=1)
            coeff=np.linalg.lstsq(design,f.ravel(),rcond=None)[0]
            residual=np.linalg.norm(design@coeff-f.ravel())
            rows.append(dict(N=n,p=p,m=m,H_coefficient=float(coeff[0].real),identity_coefficient=float(coeff[1].real),
              scalar_fit_residual=float(residual),linear_norm=float(np.linalg.norm(f)),
              commutator_norm=float(np.linalg.norm(comm[np.ix_(ids,ids)]))))
        print(n,p,rows[-(n+1):],flush=True)
    target.write_text(json.dumps(dict(status='computed',convention='E x_a x_b=delta_ab r_a/2; standard CAR; finite sizes, no DS inference',
      rows=rows,input_sha256={'code/cubic_linear_probe.py':hashlib.sha256((ROOT/'code/cubic_linear_probe.py').read_bytes()).hexdigest()}),indent=2)+'\n')
if __name__=='__main__':main()
