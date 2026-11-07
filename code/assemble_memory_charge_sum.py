"""Check inherited cancellations and charge-sum algebra without rerunning them."""
from pathlib import Path
import hashlib,json,platform
import sympy as s
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-charge-sum'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    target=OUT/'charge-sum-assembly.json'
    if target.exists():raise FileExistsError(target)
    sources=['results/susskind-memory-central-window/central-moments.json',
             'results/susskind-memory-central-norm/symbolic-limit.json',
             'results/susskind-memory-central-raw/symbolic-descendants.json',
             'results/susskind-memory-central-raw/central-assembly.json',
             'results/susskind-memory-charge-sum/global-sector-checks.json',
             'results/susskind-memory-charge-sum/growing-window-checks.json']
    prior={n:json.loads((ROOT/n).read_text()) for n in sources}
    for name,d in prior.items():
        assert d['status']=='passed',name
        for dep,h in d.get('input_sha256',{}).items():assert sha(ROOT/dep)==h,dep
    moments=prior[sources[0]]['symbolic_rows'];assert len(moments)==6
    assert sum(row['cancelled_negative_powers'] for row in moments)==21
    assert all(row['central_charge_z_cancels'] for row in moments)
    assert prior[sources[1]]['averaged_coefficients']==['0','0','4*lam**2']
    p,N,j,m,eps,t,lam=s.symbols('p N j m eps t lam',positive=True)
    assert s.summation(p-2*j,(j,0,p-1))==p
    factor=4*(m-j)*(N-m+p-j)
    assert s.factor((N+p-2*j)**2-factor-(2*m-N-p)**2)==0
    charge=s.series(s.log(s.cos(t*eps/2))/(lam*eps**2),eps,0,3)
    assert s.simplify(charge.removeO().coeff(eps,0)+t*t/(8*lam))==0
    # Conservative N^8 squared norm with N=p^2/lambda_min, A=20/lambda_min.
    exponent=16-2*lam*(20/lam);assert exponent==-24
    data=dict(status='passed',variance_completion_of_square='(N+p-2j)^2-4(m-j)(N-m+p-j)=(2m-N-p)^2',
        variance_log_envelope='2 sum_{j=0}^{p-1}(p-2j)/N=2p/N',
        charge_log_characteristic=str(charge),limiting_charge_variance='1/(4 lambda)',
        inherited_negative_moment_cancellations=21,inherited_A_norm_cancellations=2,
        outer_tail_p_exponent=int(exponent),cutoff='R_N=sqrt((20/lambda_min) log p)',
        scope='Algebra, source provenance, and applicability of polynomial identities. The logarithmic-window remainder audit and full-charge convergence proof are in manuscript, not certified by numerical scans.',
        environment=dict(python=platform.python_version(),sympy=s.__version__),
        input_sha256={n:sha(ROOT/n) for n in ['code/assemble_memory_charge_sum.py']+sources})
    target.write_text(json.dumps(data,indent=2)+'\n')
    print('charge-sum algebra and six input certificates passed')

if __name__=='__main__':main()
