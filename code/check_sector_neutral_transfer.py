"""Fresh interface checks: finite composite filtering and selected descendants."""
from datetime import datetime, timezone
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import platform
import sys
import time
import numpy as np
import sector_neutral_transfer as nt
import memory_wick_descendant_chords as dc

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/susskind-sector-neutral-transfer/interface-checks.json'


def main():
    if OUT.exists():
        raise FileExistsError(OUT)
    start=time.time()
    finite=[]
    for args in ((3,1,2,5),(4,1,2,4),(4,2,2,3)):
        finite.append(nt.filter_rows(*args))
        print(json.dumps({'finite':args,'status':'passed'}),flush=True)
    eta,V,_=dc.descendants(5)
    exact=[]
    selected=[]
    for gamma,lam in ((F(1,2),F(1,3)),(F(2,3),F(1,4))):
        q=gamma**4;s=gamma**2
        regressions,_,_=nt.regression_fields(5,gamma,lam)
        centered=[dc.fadd(field,reg,-1,f'N{j}') for j,(field,reg) in enumerate(zip(V,regressions))]
        pairs=[(1,1,0,0),(1,2,1,0),(2,2,0,0),(1,3,0,0),
               (2,3,1,0),(3,3,0,0),(1,4,1,0),(1,5,0,0)]
        for j,k,a,b in pairs:
            fields=[V[j],dc.Field(a,{dc.Term(a):dc.mono()}),
                    V[k],dc.Field(b,{dc.Term(b):dc.mono()})]
            observed=nt.field_word(fields,gamma,lam,dagger_first=True)
            old,_=dc.moment(V[j],V[k],a,b,max_slots=6)
            expected=nt.peval(old,gamma,lam)
            assert observed==expected,(gamma,lam,j,k,a,b)
            exact.append(dict(kind='two_energy_trace',gamma=str(gamma),lam=str(lam),
                              j=j,k=k,a=a,b=b,value=str(observed)))
        for j in range(1,6):
            for energy in range(0,7-j):
                energyfield=dc.Field(energy,{dc.Term(energy):dc.mono()})
                value=nt.field_word([centered[j],energyfield],gamma,lam)
                assert value==0,('regression',j,energy,value)
                exact.append(dict(kind='regression',gamma=str(gamma),lam=str(lam),
                                  j=j,energy=energy,value=str(value)))
        for j,k,a,b in pairs:
            fields=[centered[j],dc.Field(a,{dc.Term(a):dc.mono()}),
                    centered[k],dc.Field(b,{dc.Term(b):dc.mono()})]
            trace=nt.field_word(fields,gamma,lam,dagger_first=True)
            chosen=nt.field_word(fields,gamma,lam,'selected',True)
            chi=nt.field_word(fields,gamma,lam,'chi',True)
            assert chosen==(1-s*(1+lam))*trace+s*chi
            selected.append(dict(gamma=str(gamma),lam=str(lam),j=j,k=k,a=a,b=b,
                                 trace=str(trace),selected=str(chosen),chi=str(chi)))
        slope=nt.field_word([centered[1],centered[1]],gamma,lam,'selected')
        target=16*lam**3*q**2*(1+lam*s*(1+q))
        assert slope==target
        h=dc.Field(1,{dc.Term(1):dc.mono()})
        fake=(1-s)*nt.field_word([centered[1],centered[1]],gamma,lam)
        fake+=s*nt.field_word([h,h,centered[1],centered[1]],gamma,lam)
        assert slope-fake==16*lam**4*q**2*s*(1+q)>0
        selected.append(dict(kind='D_only_density_countercheck',gamma=str(gamma),lam=str(lam),
                             actual=str(slope),D_only=str(fake),difference=str(slope-fake)))
        print(json.dumps({'selected_parameter':[str(gamma),str(lam)],'status':'passed'}),flush=True)
    inputs={}
    for module in list(sys.modules.values()):
        path=getattr(module,'__file__',None)
        if path:
            path=Path(path).resolve()
            if path.parent==ROOT/'code' and path.suffix=='.py':
                inputs[path.relative_to(ROOT).as_posix()]=hashlib.sha256(path.read_bytes()).hexdigest()
    result=dict(status='passed',created_at_utc=datetime.now(timezone.utc).isoformat(),
                finite=finite,exact=exact,selected=selected,
                counts=dict(finite_rows=sum(len(x['rows']) for x in finite),
                            nonzero_filter_rows=sum(x['eta_filter_norm_squared']>1e-20 for s in finite for x in s['rows']),
                            exact_comparisons=len(exact),selected_rows=len(selected)),
                parameters=dict(C=1,formal_parameters='q=gamma^4; gamma and lambda are independent rational variables, not points on the physical curve'),
                scope='Finite full-coordinate Hermite filter/Casimir consistency and exact local-source interface checks. No finite-N pinching simulation or numerical proof of the all-degree limit.',
                environment=dict(python=platform.python_version(),numpy=np.__version__,rational='fractions.Fraction'),
                input_sha256=inputs,elapsed_seconds=time.time()-start)
    OUT.parent.mkdir(parents=True,exist_ok=True)
    with OUT.open('x',encoding='utf-8') as f:
        json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({'status':'passed','counts':result['counts'],'elapsed_seconds':result['elapsed_seconds'],
                      'sha256':hashlib.sha256(OUT.read_bytes()).hexdigest()}))


if __name__=='__main__':
    main()
