"""Fresh weight identities and reuse of immutable full-kernel quadrature nodes.

Numerical diagnostics only: no interval certification or full-spectrum limit.
Never overwrites outputs. Parameters and complete complex values are recorded.
"""
import hashlib,json,platform
from pathlib import Path
import mpmath as mp
import sympy as sp
from ds3_physical_local import (amplitude,raw_amplitude,weight,limiting_weight,
    limit_r_derivatives,labels,energy)


def enc(x):
    if isinstance(x,dict):return {k:enc(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [enc(v) for v in x]
    if isinstance(x,mp.mpc):return [mp.nstr(x.real,78),mp.nstr(x.imag,78)]
    if isinstance(x,mp.mpf):return mp.nstr(x,78)
    return x


def main():
    root=Path(__file__).resolve().parents[1]
    base=root/'results/stage4/ds3-physical-local';base.mkdir(parents=True,exist_ok=True)
    dest=base/'diagnostics.json'
    if dest.exists():raise RuntimeError('Preserve existing diagnostics')
    mp.mp.dps=85
    checks=[];points=[];derivatives=[];integrals=[]
    def check(name,value,tolerance):
        checks.append(dict(name=name,error=value,tolerance=tolerance,passed=bool(value<=tolerance)))
    jac=sp.Matrix([[1,1,1,0],[1,1,0,0],[1,0,0,0],[1,0,1,1]])
    check('absolute four-label Jacobian',abs(abs(int(jac.det()))-1),0)
    for M in [4,16,64,256]:
        for us,vs,ps,qs,rs in [('.1','-.1','.4','.4','0'),('.1','.15','.4','.5','.2'),('-.2','.1','.3','.4','-.2')]:
            u,v,p,q,r=map(mp.mpf,[us,vs,ps,qs,rs]);a,b,c,d=labels(M,u,v,r)
            val=amplitude(M,u,v,p,q,r);raw=raw_amplitude(M,u,v,p,q,r)
            rel=abs(val-raw)/abs(raw);check(f'E/W/rho M{M} {us}/{vs}/{rs}',rel,mp.mpf('1e-65'))
            points.append(dict(M=M,u=u,v=v,p=p,q=q,r=r,amplitude=val,raw=raw,relative_difference=rel,
                               energy_measure_jacobian=mp.fprod(energy(x)/x for x in [a,b,d])))
            for kind in ['crossed','exchange']:
                t=mp.mpf('1.25');analytic=limit_r_derivatives(u,v,p,q,r,t,kind)
                for j in [0,1,2]:
                    num=mp.diff(lambda z:limiting_weight(u,v,p,q,z,t,kind),r,j)
                    check(f'limiting derivative {M}/{us}/{kind}/{j}',abs(num-analytic[j]),mp.mpf('1e-70'))
                    actual=mp.diff(lambda z:weight(M,u,v,p,q,z,t,kind),r,j)
                    remainder=actual-(1+(u+v)/M)*analytic[j]
                    derivatives.append(dict(M=M,u=u,v=v,p=p,q=q,r=r,t=t,kind=kind,order=j,
                                            actual=actual,limit=analytic[j],second_scaled_residual=M*M*remainder))
    old=root/'results/stage4/ds3-uniform-local/integral-nodes.jsonl'
    nodes=[json.loads(line) for line in old.read_text().splitlines()]
    u,v,p,q,R=map(mp.mpf,['.1','-.1','.4','.4','.25'])
    for M in [64,256]:
        for order in [24,48]:
            subset=[n for n in nodes if n['M']==M and n['order']==order]
            assert len(subset)==order
            for t in map(mp.mpf,['0','1','2']):
                total=mp.mpc(0);tail=mp.mpf(0);assembled=mp.mpc(0)
                for n in subset:
                    r,w,f=map(mp.mpf,[n['r'],n['weight'],n['test']]);U=mp.mpc(*map(mp.mpf,n['value']))
                    physical=weight(M,u,v,p,q,r,t)
                    raw=raw_amplitude(M,u,v,p,q,r)*mp.exp(1j*t*(energy(M+u+v)+energy(M)-energy(M+u)-energy(M+v+r)))
                    total+=w*f*physical*U;assembled+=w*f*raw*U;tail+=abs(w*f*physical)*mp.mpf(n['tail'])
                check(f'weighted full-node assembly M{M} n{order} t{t}',abs(total-assembled),mp.mpf('1e-65'))
                target=limiting_weight(u,v,p,q,mp.mpf(0),t)
                integrals.append(dict(M=M,order=order,t=t,value=total,raw_assembly=assembled,target=target,
                                      ratio_to_target=total/target,propagated_series_tail=tail))
    files=['code/ds3_physical_local.py','code/check_ds3_physical_local.py','code/ds3_two_matter.py',
           'code/ds3_group_pairing.py','results/stage4/ds3-uniform-local/integral-nodes.jsonl']
    result=dict(status='passed' if all(c['passed'] for c in checks) else 'failed',checks=checks,
                jacobian_determinant=int(jac.det()),points=points,derivatives=derivatives,integrals=integrals,
                precision_dps=85,dependencies=dict(python=platform.python_version(),mpmath=mp.__version__,sympy=sp.__version__),
                model=dict(G=0,ell=1,threshold='fixed s>0; numerical label points all exceed .2',R='.25',
                           state='conditional initial representation label c=M; no averaging over a fixed packet'),
                scope='Local physical weights and reused full-kernel nodes, not fresh full-kernel evaluations or noncompact integrals. Derivative residuals are observations, not proof. Series tails exclude quadrature and rounding.',
                hashes={f:hashlib.sha256((root/f).read_bytes()).hexdigest() for f in files})
    dest.write_text(json.dumps(enc(result),indent=2),encoding='utf-8')
    print(json.dumps(dict(status=result['status'],checks=len(checks),points=len(points),derivatives=len(derivatives),integrals=len(integrals))))
    if result['status']!='passed':raise SystemExit(1)


if __name__=='__main__':main()
