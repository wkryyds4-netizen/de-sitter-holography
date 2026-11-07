"""Exact positive spectral matrix and the explicitly projected Taylor polynomial."""
from pathlib import Path
import hashlib,json,platform
import sympy as s

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-mixed-spectrum'


def main():
    target=OUT/'matrix-assembly.json'
    if target.exists():raise FileExistsError(target)
    W,VF,VR,V4,V2=s.symbols('W VF VR V4 V2',positive=True)
    K1,K2,K3,K4=s.symbols('K1 K2 K3 K4',positive=True)
    a,e,ss,uu=s.symbols('a e s u',real=True)
    transfer=s.Matrix([[1,0,0,0,0],[a*e,1,1,0,0],[0,0,0,1,1]])
    diag=s.diag(W*K2,VF*K1,VR*K3,V4*K4,V2*K2)
    matrix=transfer*diag*transfer.T
    expected=s.Matrix([[W*K2,a*e*W*K2,0],
                       [a*e*W*K2,VF*K1+VR*K3+a*a*e*e*W*K2,0],
                       [0,0,V4*K4+V2*K2]])
    assert matrix==expected
    det=W*K2*(VF*K1+VR*K3)*(V4*K4+V2*K2)
    assert s.factor(matrix.det()-det)==0
    left=s.Matrix([[1,-s.I*ss,-ss**2/2]])
    right=s.Matrix([1,s.I*uu,-uu**2/2])
    joint=(left*matrix*right)[0]
    formula=W*K2*(1+s.I*a*e*(uu-ss)+a*a*e*e*ss*uu)+ss*uu*(VF*K1+VR*K3)+ss**2*uu**2*(V4*K4+V2*K2)/4
    assert s.expand(joint-formula)==0
    assert s.simplify(s.conjugate(joint.xreplace({ss:uu,uu:ss}))-joint)==0
    co={(i,j):s.expand(joint).coeff(ss,i).coeff(uu,j) for i in range(3) for j in range(3)}
    assert all(co[i,j]==0 for i,j in [(0,2),(2,0),(1,2),(2,1)])
    assert s.factor(co[1,1]-matrix[1,1])==0
    assert s.factor(co[2,2]-matrix[2,2]/4)==0
    # Energy reflection (-E',-E) sends e to -e at fixed frequency.
    even=s.expand((joint+joint.subs(e,-e))/2)
    assert s.expand(even-(W*K2+ss*uu*matrix[1,1]+ss**2*uu**2*matrix[2,2]/4))==0
    names=['code/memory_mixed_assembly.py',
           'results/susskind-memory-mixed-spectrum/cubic-force-cross.json',
           'results/susskind-memory-mixed-spectrum/charged-lift.json',
           'results/susskind-memory-oriented-kernel/local-evidence.json',
           'results/susskind-memory-fluctuations/local-evidence.json']
    data=dict(status='passed',matrix=[[str(x) for x in row] for row in matrix.tolist()],
              positive_factorization=dict(transfer=transfer.tolist(),diagonal=list(map(str,diag.diagonal()))),
              determinant=str(det),projected_polynomial_energy_density=str(s.expand(joint)),
              coefficients={f'{i},{j}':str(v) for (i,j),v in co.items()},
              ordinary_time_kernel='W G_q2(T)+s*u C_L(T)+s^2*u^2 C_X(T)/4',
              scope='Phi(s)=Y+i*s*L-s^2*X/2 uses the separately projected second coefficient X. It is not the raw memory Taylor polynomial and has no finite-time remainder claim.',
              environment=dict(python=platform.python_version(),sympy=s.__version__),
              input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names})
    # SymPy integers in the transfer matrix are made explicit JSON integers/strings.
    data['positive_factorization']['transfer']=[[str(v) for v in row] for row in transfer.tolist()]
    target.write_text(json.dumps(data,indent=2)+'\n')
    print('passed positive 3x3 factorization, determinant, and projected polynomial coefficients',flush=True)


if __name__=='__main__':main()
