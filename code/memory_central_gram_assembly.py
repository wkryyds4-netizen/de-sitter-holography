"""Symbolic central Gram projection after the polarization proof."""
from pathlib import Path
import hashlib,json,platform
import sympy as s

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-central-gram'

def main():
    target=OUT/'gram-assembly.json'
    if target.exists():raise FileExistsError(target)
    lam,q,Js=s.symbols('lambda q J_s',positive=True)
    z=1+q;b=2*lam*(1-q);l=2*lam*(1+q)
    ss=4*lam**2*(1+q)+2*lam*(1-q);aa=4*lam+4*lam**2*(1+q)
    v=s.factor(ss-b*b/z);delta=s.factor(aa-ss)
    G=s.Matrix([[z,0,b,b],[0,l,0,0],[b,0,ss,ss],[b,0,ss,aa]])
    assert s.factor(G.det()-z*l*v*delta)==0
    bo=v/(v+delta);bb=b*delta/(z*(v+delta))
    residual=s.Matrix([-bb,0,1,-bo])
    for j in [0,1,3]:assert s.factor((G*residual)[j])==0
    variance=s.factor((residual.T*G*residual)[0])
    assert s.factor(variance-v*delta/(v+delta))==0
    core=s.factor(64*lam**2*q**4*variance)
    physical=s.factor(core*Js**8/lam**4)
    small=s.limit(physical.subs(q,s.exp(-4*lam)),lam,0,dir='+')
    assert small==1024*Js**8
    data=dict(status='passed',basis=['B','L_input','S','O'],gram_over_C_squared=[[str(t) for t in row] for row in G.tolist()],
        variance_S_perp=str(v),variance_O_minus_S=str(delta),projection_O=str(s.factor(bo)),projection_B=str(s.factor(bb)),
        residual_variance_over_C_squared=str(variance),core_residual_variance_over_C_fourth=str(core),
        physical_core_variance=str(physical),ordered_small_lambda_limit=str(small),
        positivity_conditions='lambda>0 and q=exp(-4 lambda) in (0,1); v and Delta strictly positive.',
        scope='Symbolic assembly of the proven equal-time central-window Gram. Core uses the previously proven strong W2 expansion. Does not extend energy-resolved spectra, integrate all charges, or prove the optional exact charged Casimir formula.',
        environment=dict(python=platform.python_version(),sympy=s.__version__),
        input_sha256={'code/memory_central_gram_assembly.py':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    target.write_text(json.dumps(data,indent=2)+'\n');print('central four Gram and positive core residual passed; small lambda:',small)

if __name__=='__main__':main()
