"""Fresh exact joint-X regression, residual Gram and mixed spectral checks.

All arithmetic is rational, at independent formal q and lambda. C=1.
No predecessor main is run, and no existing scientific output is overwritten.
"""
from collections import Counter
from datetime import datetime, timezone
from fractions import Fraction as F
from itertools import combinations_with_replacement
from pathlib import Path
import hashlib
import json
import platform
import time

import numpy as np
import memory_wick_descendant_chords as neutral_source
from memory_joint_regression import JointModel

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-joint-regression/joint-regression-checks.json'
PARAMETERS={
    'models':[{'q':'1/3','lambda':'1/4','maximum_G_degree':5,'maximum_V_degree':4},
              {'q':'1/2','lambda':'1/2','maximum_G_degree':3,'maximum_V_degree':3}],
    'C':1,'maximum_raw_slots':10,'charge_test_powers':[0,1,2,3],
    'Gram_energy_insertions':[[0,0],[1,0],[0,1],[1,1],[2,0],[0,2]],
    'no_bridge_energy_insertions':[[0,0],[1,0],[0,1],[1,1],[2,0],[0,2]],
    'mixed_energy_insertions':[[0,0],[1,0],[0,1],[1,1]],
    'formal_parameter_scope':'q and lambda are independent rational variables; these points are not asserted to satisfy q=exp(-4 lambda)',
}
PARENTS=['results/susskind-memory-oriented-transfer/body-geometry-checks.json',
         'results/susskind-memory-wick-descendants/descendant-checks.json',
         'results/susskind-memory-wick-descendants/source-regression-checks.json']


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':')).encode()).hexdigest()


def evaluate(poly,m):
    return sum((v*m.lam**i*m.q**j for (i,j),v in poly.items()),F(0))


def exact_pair(m,left,right,a=0,b=0,no_bridge=False):
    slots=max((t.length for t in left),default=0)+max((t.length for t in right),default=0)+a+b
    assert slots<=PARAMETERS['maximum_raw_slots']
    return m.moment(left,right,a,b,no_bridge)


def checks(spec):
    m=JointModel(spec['q'],spec['lambda']);nmax=spec['maximum_G_degree'];vmax=spec['maximum_V_degree']
    _,old_V,_=neutral_source.descendants(vmax)
    V={j:m.neutral(old_V[j]) for j in range(1,vmax+1)}
    R={n:m.residual(n) for n in range(2,nmax+1)}
    S={n:m.resummed(n) for n in range(2,nmax+1)}
    counts=Counter();projection=[];grams=[];b0=[];mixed=[];calibration=[];neutral=[]
    for n in range(2,nmax+1):
        for k in range(6):
            if n+k>8 or (n+k)%2:continue
            for r in PARAMETERS['charge_test_powers']:
                test=m.charge(m.H(k),r)
                raw=exact_pair(m,test,m.G(n));reg=exact_pair(m,test,m.g(n))
                assert raw==reg,(spec,'projection',n,k,r,raw,reg)
                projection.append({'n':n,'energy_power':k,'charge_power':r,
                                   'raw_joint_moment':str(raw),'regression_joint_moment':str(reg),'exact':True})
                counts['joint_projection_tests']+=1
                counts['nonzero_joint_projection_tests']+=raw!=0
    for n,k in combinations_with_replacement(range(2,nmax+1),2):
        for a,b in PARAMETERS['Gram_energy_insertions']:
            if n+k+a+b>PARAMETERS['maximum_raw_slots'] or (n+k+a+b)%2:continue
            rr=exact_pair(m,R[n],R[k],a,b)
            rs=exact_pair(m,R[n],S[k],a,b)
            sr=exact_pair(m,S[n],R[k],a,b)
            ss=exact_pair(m,S[n],S[k],a,b)
            assert rr==rs==sr==ss,(spec,'Gram',n,k,a,b,rr,rs,sr,ss)
            grams.append({'n':n,'m':k,'left_energy':a,'right_energy':b,
                          'RR':str(rr),'RS':str(rs),'SR':str(sr),'SS':str(ss),'exact':True})
            counts['residual_Gram_layouts']+=1;counts['residual_Gram_equalities']+=3
            counts['nonzero_residual_Gram_layouts']+=rr!=0
    for n,k in combinations_with_replacement(range(2,min(nmax,4)+1),2):
        for a,b in PARAMETERS['no_bridge_energy_insertions']:
            if n+k+a+b>PARAMETERS['maximum_raw_slots'] or (n+k+a+b)%2:continue
            raw=exact_pair(m,m.G(n),m.G(k),a,b,True)
            left=exact_pair(m,m.g(n),m.G(k),a,b,True)
            right=exact_pair(m,m.G(n),m.g(k),a,b,True)
            both=exact_pair(m,m.g(n),m.g(k),a,b,True)
            assert raw-left-right+both==0,(spec,'b0',n,k,a,b)
            b0.append({'n':n,'m':k,'left_energy':a,'right_energy':b,
                       'GG':str(raw),'gG':str(left),'Gg':str(right),'gg':str(both),
                       'centered_group_sum':'0','exact':True})
            counts['no_bridge_layouts']+=1;counts['nonzero_raw_no_bridge_layouts']+=raw!=0
    for n in range(2,nmax+1):
        for j in range(1,vmax+1):
            for a,b in PARAMETERS['mixed_energy_insertions']:
                if n+j+a+b>PARAMETERS['maximum_raw_slots'] or (n+j+a+b)%2:continue
                raw=exact_pair(m,R[n],V[j],a,b);summed=exact_pair(m,S[n],V[j],a,b)
                assert raw==summed,(spec,'mixed',n,j,a,b,raw,summed)
                mixed.append({'n':n,'neutral_V_degree':j,'left_energy':a,'right_energy':b,
                              'whole_residual_V':str(raw),'resummed_J_V':str(summed),'exact':True})
                counts['mixed_neutral_layouts']+=1;counts['nonzero_mixed_neutral_layouts']+=raw!=0
    # An independent sealed neutral engine checks mark/OU domain conversion.
    for i,j in combinations_with_replacement(range(1,min(vmax,3)+1),2):
        for a,b in ((0,0),(1,1)):
            if (i+j+a+b)%2:continue
            prior,_=neutral_source.moment(old_V[i],old_V[j],a,b,max_slots=8)
            old=evaluate(prior,m);new=exact_pair(m,V[i],V[j],a,b)
            assert old==new,(spec,'neutral conversion',i,j,a,b)
            neutral.append({'i':i,'j':j,'left_energy':a,'right_energy':b,
                            'sealed_algorithm':str(old),'new_algorithm':str(new),'exact':True})
            counts['neutral_algorithm_crosschecks']+=1
    # A0=J2; B2=2C(V2-{H,V1}) after static removal, here C=1.
    # A0 is already orthogonal to the static subtraction of B2.
    A0=m.J(2);B1=m.scale(V[1],2)
    B2=m.scale(m.add(V[2],m.scale(m.hleft(V[1]),-1),m.scale(m.hright(V[1]),-1)),2)
    w=4*m.lam+16*m.lam*m.lam*m.q/(1+m.q)
    eta=-4*m.lam*m.q+16*m.lam*m.lam*m.q/(1+m.q)
    prefactor=4*m.lam*m.q*m.q*(w+eta)
    for degree in (0,2,4):
        for a in range(degree+1):
            b=degree-a
            direct=exact_pair(m,A0,B2,a,b)
            mu=m.mehler_moment(a,b,m.q*m.q)
            reference=prefactor*mu
            assert direct==reference and reference>0,(spec,'A0B2',a,b,direct,reference)
            calibration.append({'channel':'A0B2','left_energy':a,'right_energy':b,
                                'direct':str(direct),'Mehler_moment_q_squared':str(mu),
                                'known_prefactor':str(prefactor),'reference':str(reference),'exact':True})
            counts['nonzero_A0B2_Mehler_checks']+=1
    for degree in (1,3):
        for a in range(degree+1):
            b=degree-a
            value=exact_pair(m,A0,B1,a,b)
            assert value==0
            calibration.append({'channel':'A0B1','left_energy':a,'right_energy':b,'direct':'0','exact':True})
            counts['zero_A0B1_checks']+=1
    energy_g2={term:value for term,value in m.g(2).items() if term.xpower==0}
    energy_residual=m.add(m.G(2),m.scale(energy_g2,-1))
    energy_b0=exact_pair(m,energy_residual,energy_residual,no_bridge=True)
    joint_b0=exact_pair(m,R[2],R[2],no_bridge=True)
    full_norm=exact_pair(m,R[2],R[2])
    assert energy_b0==4*m.lam and joint_b0==0 and full_norm==w
    counts.update({'charge_static_diagnostic':1,'G2_joint_norm_calibration':1})
    diagnostic={'energy_only_no_bridge':'retained nonzero charge background',
                'energy_only_b0_value':str(energy_b0),'known_charge_variance':str(4*m.lam),
                'joint_centered_b0_value':str(joint_b0),'joint_G2_norm':str(full_norm),'old_w':str(w)}
    return {'parameters':spec,'counts':dict(counts),'diagram_statistics':dict(m.stats),
            'candidate_fields':[{'n':n,'g_n':m.serialize(m.g(n)),'resummed_J':m.serialize(S[n])}
                                for n in range(2,nmax+1)],
            'joint_projection':projection,'residual_Grams':grams,'no_bridge':b0,
            'mixed_neutral':mixed,'neutral_algorithm_crosschecks':neutral,
            'low_order_spectral_calibration':calibration,'charge_static_diagnostic':diagnostic}


def main():
    if OUT.exists():raise FileExistsError(OUT)
    started=time.perf_counter();parents={}
    for name in PARENTS:
        data=json.loads((ROOT/name).read_text(encoding='utf-8'));assert data['status']=='passed'
        for input_name,want in data['input_sha256'].items():assert sha(ROOT/input_name)==want,(name,input_name)
        parents[name]={'sha256':sha(ROOT/name),'input_hashes_verified':True}
    rows=[];counts=Counter();statistics=Counter()
    for spec in PARAMETERS['models']:
        row=checks(spec);rows.append(row);counts.update(row['counts']);statistics.update(row['diagram_statistics'])
        print(json.dumps({'completed_model':spec,'counts':row['counts']},sort_keys=True),flush=True)
    inputs=['code/memory_joint_regression.py','code/check_memory_joint_regression.py',
            'code/memory_wick_descendant_chords.py','code/memory_body_source_chords.py',
            'code/memory_oriented_marked_chords.py','code/memory_size_marked_chords.py',
            'code/singlet_coefficient.py','code/singlet_multimark.py','code/complex_syk_dictionary.py']+PARENTS
    result={'status':'passed','schema_version':1,'created_utc':datetime.now(timezone.utc).isoformat(),
            'scope':'Exact fixed-degree formal graph checks of candidate joint charge/energy regression and resummed residual Grams, including mixed neutral and nonzero A0B2. Uses established joint-X graph law and known Mehler moments; finite tests do not replace all-degree analytic, uniform-integrability or pinching proofs.',
            'normalization':'C=1; X centered Gaussian with variance 1/(4 lambda). q and lambda are independent formal rational parameters.',
            'parameters':PARAMETERS,'parameter_sha256':digest(PARAMETERS),
            'environment':{'python':platform.python_version(),'numpy':np.__version__,'platform':platform.platform()},
            'input_sha256':{name:sha(ROOT/name) for name in inputs},'parents_verified':parents,
            'counts':dict(counts),'diagram_statistics':dict(statistics),'models':rows,
            'Mehler_source':'Existing q-Hermite/Mehler density: moment(E^a Eprime^b)=sum_d a_ad a_bd r^d [d]_q!; see literature/susskind-memory-oriented-kernel-sources.md and notes/susskind-memory-two-component.md.',
            'elapsed_seconds':time.perf_counter()-started}
    OUT.parent.mkdir(parents=True,exist_ok=True)
    with OUT.open('x',encoding='utf-8') as stream:
        json.dump(result,stream,ensure_ascii=False,indent=2);stream.write('\n')
    print(json.dumps({'status':'passed','counts':dict(counts),'diagram_statistics':dict(statistics),
                      'output_sha256':sha(OUT),'parameter_sha256':result['parameter_sha256'],
                      'elapsed_seconds':result['elapsed_seconds']},sort_keys=True))


if __name__=='__main__':main()
