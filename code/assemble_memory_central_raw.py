"""Assemble central conserved shift and spectral matrix; reuse old certificates."""
from pathlib import Path
import hashlib,json,platform
import sympy as s

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-central-raw'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    target=OUT/'central-assembly.json'
    if target.exists():raise FileExistsError(target)
    names=['results/susskind-memory-central-raw/symbolic-descendants.json',
           'results/susskind-memory-central-raw/central-filter-checks.json',
           'results/susskind-memory-raw-spectrum/symbolic-spectrum.json',
           'results/susskind-memory-mixed-spectrum/cubic-force-cross.json',
           'results/susskind-memory-mixed-spectrum/charged-lift.json',
           'results/susskind-filtered-spectrum/charged-bridge.json',
           'results/susskind-memory-central-spectrum/spectrum-assembly.json']
    data={n:json.loads((ROOT/n).read_text()) for n in names}
    for name,d in data.items():
        assert d['status']=='passed',name
        for dep,h in d.get('input_sha256',{}).items():assert sha(ROOT/dep)==h,dep
    q,lam,x,C=s.symbols('q lam x C');old=data[names[2]]
    co={k:s.sympify(v) for k,v in old['coefficients'].items()}
    shift=s.sympify(data[names[0]]['raw_Q4_B_shift'])
    assert s.factor(shift-2*lam*x*co['B_leading'])==0
    # Four-field positive decomposition for (Y, L, U2).
    w,f,r,z=s.symbols('w f r z',positive=True)
    e,H,a,bF,bR=s.symbols('e H a bF bR',real=True)
    M=s.Matrix([[1,0,0,0],[a*e,1,1,0],[H,bF*e,bR*e,1]])
    mat=M*s.diag(w,f,r,z)*M.T
    determinant=w*((f+r)*z+e*e*f*r*(bF-bR)**2)
    assert s.factor(mat.det()-determinant)==0
    assert mat[0,2]==w*H
    assert s.expand(mat[1,2]-(a*e*w*H+bF*e*f+bR*e*r))==0
    # Exact forward/reverse lift coefficients have the correct central sizes.
    eps=s.symbols('eps',positive=True);pp=1/eps;mm=1/(2*lam*eps**2)+x/eps;d=mm-pp+1
    limits=[s.limit(pp/d,eps,0),s.limit(1/(eps*d),eps,0),s.limit(pp/(eps*d),eps,0)]
    assert limits==[0,0,2*lam]
    result=dict(status='passed',central_counterterm_extra_B_coefficient=str(C*shift),
        identical_dynamic_coefficients={k:str(co[k]) for k in ['bR','bF','bJ','bD','kappa']},
        positive_transfer=[[str(v) for v in row] for row in M.tolist()],
        determinant=str(determinant),mixed_Y_U2=str(mat[0,2]),mixed_L_U2=str(mat[1,2]),
        scope='Algebraic assembly and input provenance, with central convergence proved in manuscript. No fresh energy-kernel scan; no finite-memory remainder claim.',
        reused_certificates=names[2:],environment=dict(python=platform.python_version(),sympy=s.__version__),
        input_sha256={n:sha(ROOT/n) for n in ['code/assemble_memory_central_raw.py']+names})
    target.write_text(json.dumps(result,indent=2)+'\n')
    print('central counterterm and positive matrix passed; seven prior/current certificates verified')

if __name__=='__main__':main()
