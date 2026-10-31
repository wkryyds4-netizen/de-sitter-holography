"""Local test-function integrals: full kernel at finite M and sine asymptotic separately."""
import hashlib,json,platform,time
from pathlib import Path
import mpmath as mp
from ds3_uniform_local import evaluate
from ds3_transverse_phase import scattering_phase
from run_ds3_collision import clean


def main():
    root=Path(__file__).resolve().parents[1];base=root/'results/stage4/ds3-uniform-local'
    dest=base/'integrals.json';raw=base/'integral-nodes.jsonl'
    if dest.exists() or raw.exists():raise RuntimeError('Preserve prior integrals')
    mp.mp.dps=60;R=mp.mpf('.25');full=[];asym=[]
    for M in [64,256]:
        for order in [24,48]:
            nodes,weights=mp.gauss_quadrature(order,'legendre');sums=[mp.mpc(0)]*3;tails=[mp.mpf(0)]*3
            for j in range(order):
                r=R*nodes[j];w=R*weights[j];f=(1-(r/R)**2)**3
                row=evaluate(M,'.1','-.1','.4','.4',str(r),R=str(R))
                for k in range(3):sums[k]+=w*f*r**k*row['value'];tails[k]+=w*f*abs(r)**k*row['tail']
                node=clean(dict(M=M,order=order,index=j,r=r,weight=w,test=f,value=row['value'],tail=row['tail']))
                with raw.open('a',encoding='utf-8') as out:out.write(json.dumps(node)+'\n');out.flush()
            full.append(clean(dict(M=M,order=order,integrals=sums,propagated_series_tails=tails)))
            print('full',M,order,[mp.nstr(x,10) for x in sums],flush=True)
    for L in [4,8,16,32,64,128]:
        for order in [64,128]:
            nodes,weights=mp.gauss_quadrature(order,'legendre');sums=[mp.mpf(0)]*3
            for j in range(order):
                r=R*nodes[j];w=R*weights[j];f=(1-(r/R)**2)**3
                K=mp.sin(L*r-scattering_phase('.1','-.1','.4','.4',r))/(mp.pi*r)
                for k in range(3):sums[k]+=w*f*r**k*K
            asym.append(clean(dict(logM=L,order=order,integrals=sums)))
        print('sine',L,[mp.nstr(x,10) for x in sums],flush=True)
    paths=['code/integrate_ds3_uniform_local.py','code/ds3_uniform_local.py','code/ds3_transverse_phase.py']
    result=dict(status='completed',R=str(R),tests=['(1-(r/R)^2)^3','r*(1-(r/R)^2)^3','r^2*(1-(r/R)^2)^3'],
                targets=['1','0','0'],full_kernel=full,sine_model=asym,dps=60,
                scope='Compact dr tests only. Gauss refinement is empirical; propagated series tails exclude quadrature and rounding. Large logM rows evaluate only the proved sine asymptotic.',
                hashes={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in paths})
    dest.write_text(json.dumps(result,indent=2),encoding='utf-8')


if __name__=='__main__':main()
