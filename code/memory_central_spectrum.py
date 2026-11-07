"""Central oriented-kernel assembly, with explicit reused certificate hashes."""
from pathlib import Path
import hashlib,json,platform
import sympy as s

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-central-spectrum'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    target=OUT/'spectrum-assembly.json'
    if target.exists():raise FileExistsError(target)
    sources=['results/susskind-memory-oriented-kernel/reverse-cross.json',
             'results/susskind-memory-oriented-kernel/oriented-marks.json',
             'results/susskind-memory-oriented-kernel/spectral-checks.json',
             'results/susskind-memory-quartic/cross-kernel.json',
             'results/susskind-memory-quartic/quartic.json']
    prior={name:json.loads((ROOT/name).read_text()) for name in sources}
    for name,data in prior.items():
        assert data['status']=='passed',name
        for dep,h in data.get('input_sha256',{}).items():assert sha(ROOT/dep)==h,dep
    rev=prior[sources[0]];assert rev['certificate']['templates']==128
    assert all(r['vanishes'] for r in rev['certificate']['records'])
    assert not rev['certificate']['counterexamples']
    lam,q,C=s.symbols('lambda q C',positive=True)
    v=2*lam*(1-q)+16*lam**2*q/(1+q);d=2*lam*(1+q);rho=v/(v+d)
    w=4*C**2*lam+16*C**2*lam**2*q/(1+q)
    cross=-4*C**2*lam*q+16*C**2*lam**2*q/(1+q)
    M=s.Matrix([[w,cross],[cross,w]])
    S=s.Matrix([s.Rational(1,2),s.Rational(1,2)]);T=S-s.Matrix([rho,0])
    assert s.factor((S.T*M*S)[0]-C**2*v)==0
    assert s.factor((M*T)[0])==0
    vt=s.factor((T.T*M*T)[0]);assert s.factor(vt-C**2*v*d/(v+d))==0
    v2=s.factor((8*lam*C*q**2)**2*vt)
    old=s.sympify(rev['V2'],locals={'lam':lam,'C':C,'q':q})
    assert s.factor(v2-old)==0
    curvature=s.factor(2*C*(1-q*q))
    data=dict(status='passed',oriented_matrix=[[str(t) for t in row] for row in M.tolist()],
        spectral_parameter='q^2',T_coefficients=list(map(str,T)),T_variance=str(vt),projected_W2_core_variance=str(v2),
        normalized_curvature=str(curvature),cross_measures_zero=['T--Y_minus','Q4--Y_minus','Q4--Y_plus','Q4--T','projected_W2_core--Q4'],
        reused_reverse_templates=128,reused_reverse_words=6072,
        scope='Uses the new uniform central trace theorem and general-charge orientation bridge. Reuses existing marked graph kernels and endpoint certificates, verifying their recorded source hashes without rerunning them. Proves the projected W2 core spectrum; does not yet restore all central raw U2 descendants or filtered F crosses.',
        environment=dict(python=platform.python_version(),sympy=s.__version__),
        input_sha256={'code/memory_central_spectrum.py':sha(Path(__file__)),**{name:sha(ROOT/name) for name in sources}})
    target.write_text(json.dumps(data,indent=2)+'\n');print('central oriented spectrum and projected W2 core assembled; reused certificates unchanged')

if __name__=='__main__':main()
