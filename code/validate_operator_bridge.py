"""Finite KMS oscillator checks of the observer Weyl/field bridge.

These numerical checks do not prove continuum localization or affiliation.
Any clock probability supported strictly in (omega,2*omega) yields the same
intermediate projection used here; no delta-function clock state is needed.
"""
from pathlib import Path
import hashlib, json, platform
import numpy as np
import scipy
from scipy.linalg import eigh, solve

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/stage2'

def main():
    N=48; omega=1.
    n=np.arange(N)
    probabilities=np.exp(-omega*n);probabilities/=probabilities.sum()
    Omega=np.diag(np.sqrt(probabilities)).astype(complex)
    A=np.diag(np.sqrt(np.arange(1,N)),1)
    A=A+A.T
    variance=np.sum(probabilities*np.diag(A@A))
    A=A/np.sqrt(variance)
    energy=omega*(n[:,None]-n[None,:])
    P=(energy<=1.5*omega)
    evals,evecs=eigh(A)
    def functional(fun):
        return (evecs*fun(evals))@evecs.T
    def translated(matrix,t):
        return matrix*np.exp(1j*energy*t)
    def expectation(base,times):
        state=Omega.copy()
        for t in reversed(times):
            state=(translated(base,t)@state)*P
        return complex(np.vdot(Omega,state))
    def spectral_moment(times):
        a=1/(1+np.exp(-omega));b=1-a
        if len(times)==2:
            return a*np.exp(1j*omega*(times[1]-times[0]))+b*np.exp(-1j*omega*(times[1]-times[0]))
        t1,t2,t3,t4=times
        total=0j
        for x,px in [(omega,a),(-omega,b)]:
            for y,py in [(omega,a),(-omega,b)]:
                q=[max(0,x,y),max(x,0)+max(y,0),max(0,x+max(y,0))]
                phases=[x*(t2-t1)+y*(t4-t3),x*(t3-t1)+y*(t4-t2),x*(t4-t1)+y*(t3-t2)]
                total+=px*py*sum(np.exp(1j*v) for v,cut in zip(phases,q) if cut<=1.5*omega)
        return total
    tasks={'two':[.3,-.2],'four_equal':[0,0,0,0],'four_unequal':[.2,-.1,.5,0]}
    exact={key:expectation(A,times) for key,times in tasks.items()}
    checks=[]
    def check(name,passed,**kwargs):
        checks.append({'name':name,'passed':bool(passed),**kwargs})
    for key,times in tasks.items():
        error=abs(exact[key]-spectral_moment(times))
        check('projected Wick measure '+key,error<2e-12,absolute_error=error)
    # The finite matrix state is KMS but its truncated ladder is not an exact CCR field.
    # Double the level cutoff to check moments empirically, not by probability-tail alone.
    Ncheck=2*N
    nn=np.arange(Ncheck)
    rr=np.exp(-omega*nn);rr/=rr.sum()
    ref=np.diag(np.sqrt(rr)).astype(complex)
    large=np.diag(np.sqrt(np.arange(1,Ncheck)),1)
    large=large+large.T
    large/=np.sqrt(np.sum(rr*np.diag(large@large)))
    energies=omega*(nn[:,None]-nn[None,:]);project=(energies<=1.5*omega)
    cutoff_changes={}
    for key,times in tasks.items():
        state=ref.copy()
        for t in reversed(times):
            state=((large*np.exp(1j*energies*t))@state)*project
        change=abs(complex(np.vdot(ref,state))-exact[key])
        cutoff_changes[key]=change
        check('48 to 96 level convergence '+key,change<2e-12,absolute_change=change,
              scope='Empirical cutoff convergence for this moment; not an all-observable bound.')
    probe=(A@Omega)*P
    probe/=np.linalg.norm(probe)
    leakage_vector=(A@probe)*(~P)
    quadratic_difference=-(A@leakage_vector)*P/2
    quadratic_norm_loss=float(np.linalg.norm(leakage_vector)**2)
    rows=[]
    for eps in [.4,.2,.1,.05,.025,.0125]:
        sine=functional(lambda values:np.sin(eps*values))/eps
        unitary=functional(lambda values:np.exp(1j*eps*values))
        moment_errors={key:abs(expectation(sine,times)-exact[key]) for key,times in tasks.items()}
        reference_compressed=(unitary@Omega)*P
        reference_exponential=np.zeros_like(Omega)
        compressed_action=(unitary@probe)*P
        exponential_action=np.zeros_like(Omega)
        r_sine=np.zeros_like(Omega);r_field=np.zeros_like(Omega)
        for column in range(N):
            dim=min(N,column+2)
            block=A[:dim,:dim]
            vals,vecs=eigh(block)
            expblock=(vecs*np.exp(1j*eps*vals))@vecs.T
            exponential_action[:dim,column]=expblock@probe[:dim,column]
            reference_exponential[:dim,column]=expblock@Omega[:dim,column]
            vector=Omega[:dim,column]
            r_sine[:dim,column]=solve(sine[:dim,:dim]-1j*np.eye(dim),vector)
            r_field[:dim,column]=solve(block-1j*np.eye(dim),vector)
        separation=float(np.linalg.norm(compressed_action-exponential_action))
        leakage=float(1-np.linalg.norm(compressed_action)**2)
        resolvent_error=float(np.linalg.norm(r_sine-r_field))
        leading_relative_error=float(np.linalg.norm(compressed_action-exponential_action-eps**2*quadratic_difference)/(eps**2*np.linalg.norm(quadratic_difference)))
        rows.append({'epsilon':eps,'moment_errors':moment_errors,
                     'resolvent_vector_error':resolvent_error,
                     'compressed_Weyl_vs_exp_compressed_field':separation,
                     'compressed_Weyl_norm_loss':leakage,
                     'leading_difference_relative_residual':leading_relative_error,
                     'norm_loss_coefficient_ratio':leakage/(eps**2*quadratic_norm_loss),
                     'reference_vector_Weyl_difference':float(np.linalg.norm(reference_compressed-reference_exponential)),
                     'reference_vector_norm_loss':float(1-np.linalg.norm(reference_compressed)**2)})
        check(f'compressed Weyl differs at {eps}',separation>0,
              vector_difference=separation)
        check(f'compressed Weyl loses norm at {eps}',leakage>0,norm_loss=leakage)
    check('boundary probe second-order compression defect',rows[-1]['leading_difference_relative_residual']<.03,
          relative_residual=rows[-1]['leading_difference_relative_residual'])
    check('boundary probe norm loss coefficient',abs(rows[-1]['norm_loss_coefficient_ratio']-1)<.002,
          coefficient_ratio=rows[-1]['norm_loss_coefficient_ratio'])
    for key in tasks:
        errors=[row['moment_errors'][key] for row in rows]
        ratios=[errors[i]/errors[i+1] for i in range(len(errors)-1)]
        check('quadratic mixed-difference convergence '+key,
              all(errors[i+1]<errors[i] for i in range(len(errors)-1)) and
              all(3.8<r<4.2 for r in ratios[-3:]),ratios=ratios)
    ratios=[rows[i]['resolvent_vector_error']/rows[i+1]['resolvent_vector_error'] for i in range(len(rows)-1)]
    check('quadratic resolvent convergence on reference vector',all(3.8<r<4.2 for r in ratios[-3:]),ratios=ratios)
    code_hash=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result={'status':'passed' if all(c['passed'] for c in checks) else 'failed',
            'checks':len(checks),'failed':sum(not c['passed'] for c in checks),
            'scope':'Finite-matrix KMS GNS versus infinite-oscillator Wick reference; projected Weyl central differences, resolvents on one reference vector, and nonmultiplicativity. Truncated ladders do not obey exact CCR. These are supporting checks, not continuum proofs.',
            'cutoff_moment_changes':cutoff_changes,
            'probe':'normalized projected A acting on the thermal reference; the reference alone cancels the leading defect',
            'parameters':{'oscillator_levels':N,'omega':omega,'projection_threshold':1.5,'clock_support':'any L2 wave amplitude within (1,2)',
                          'omitted_thermal_probability':float(np.exp(-N*omega)),
                          'tail_note':'This probability tail is not a bound on unbounded field moments; separate 48-to-96 empirical moment checks are included.'},
            'versions':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},
            'code_sha256':{Path(__file__).name:code_hash},'tasks':tasks,
            'projected_moments':{k:[v.real,v.imag] for k,v in exact.items()},
            'rows':rows,'check_details':checks}
    OUT.mkdir(exist_ok=True)
    (OUT/'operator-bridge-validation.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:result[k] for k in ['status','checks','failed','projected_moments']},indent=2))
    print(json.dumps(rows[-1],indent=2))
    if result['failed']:raise SystemExit(1)

if __name__=='__main__':
    main()
