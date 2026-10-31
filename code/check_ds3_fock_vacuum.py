"""New vacuum-return quadratures and an independent finite Fock-path check."""
import hashlib,json,platform
from pathlib import Path
import mpmath as mp
import numpy as np
from ds3_fock_vacuum import raw_rate,return_block,block_bound
from ds3_two_matter import tilted
from run_ds3_collision import clean

def finite_fock():
    states=[(i,n-i) for n in range(5) for i in range(n+1)];index={v:k for k,v in enumerate(states)}
    V=[np.array([[0,1,1j],[.3j,.2,0],[.4,0,-.1j]],complex)/3,
       np.array([[.2j,0,.1],[1,0,.7j],[0,.4,0]],complex)/3]
    energies=np.sqrt(1+np.array([.4,1.1,2.3])**2);dim=len(states)*3
    a=[]
    for j in range(2):
        op=np.zeros((len(states),len(states)),complex)
        for st,k in index.items():
            target=list(st);target[j]+=1
            if tuple(target) in index:op[index[tuple(target)],k]=np.sqrt(st[j]+1)
        a.append(op)
    def creation(t):
        phase=np.exp(-1j*energies*t)
        return sum(np.kron(op,phase[:,None]*v*phase.conj()[None,:]) for op,v in zip(a,V))
    def field(t):
        c=creation(t);return c+c.conj().T
    masks={n:np.array([3*k+j for k,st in enumerate(states) if sum(st)==n for j in range(3)]) for n in range(5)}
    f=np.zeros(dim,complex);f[:3]=np.array([1,1j,2],complex)/np.sqrt(6)
    normV=float(np.linalg.eigvalsh(sum(v.conj().T@v for v in V))[-1]);bounds=[]
    C=creation(.7)
    for n in range(4):
        norm=float(np.linalg.norm(C[np.ix_(masks[n+1],masks[n])],ord=2))
        bounds.append(dict(N=n,norm=norm,bound=np.sqrt((n+1)*normV),passed=norm<=np.sqrt((n+1)*normV)+1e-14))
    rows=[]
    for ts in [(0,0,0,0),(.2,-.7,1.3,.4),(1,0,1,0),(-.5,.8,2,-1)]:
        t1,t2,t3,t4=ts;c1,c2,c3,c4=map(creation,ts)
        direct=np.vdot(f,field(t1)@field(t2)@field(t3)@field(t4)@f)
        vac=np.vdot(f,c1.conj().T@c2@c3.conj().T@c4@f)
        two=np.vdot(c2@c1@f,c3@c4@f)
        rows.append(dict(times=ts,full=[direct.real,direct.imag],vacuum=[vac.real,vac.imag],two=[two.real,two.imag],
            difference=abs(direct-vac-two),passed=abs(direct-vac-two)<1e-14))
    comm=[]
    for t in [.1,.7,1.4]:
        v=(field(t)@field(0)-field(0)@field(t))@f
        full=float(np.vdot(v,v).real);v0=float(np.vdot(v[masks[0]],v[masks[0]]).real);v2=float(np.vdot(v[masks[2]],v[masks[2]]).real)
        comm.append(dict(t=t,full=full,vacuum=v0,two=v2,passed=abs(full-v0-v2)<1e-14))
    return dict(scope='Finite observer matrices and two oscillator modes; tests Fock normalization and particle-number paths only, not de Sitter group averaging.',
        particle_cap=4,observer_energies=energies.tolist(),V=[[[[z.real,z.imag] for z in row] for row in v] for v in V],
        bounds=bounds,fourpoint=rows,commutators=comm)

def main():
    root=Path(__file__).resolve().parents[1];base=root/'results/stage4/ds3-fock-vacuum'
    base.mkdir(parents=True,exist_ok=True);dest=base/'diagnostics.json';raw=base/'vacuum-runs.jsonl'
    if dest.exists() or raw.exists():raise RuntimeError('Preserve existing evidence')
    mp.mp.dps=60;rates=[];rows=[];checks=[]
    for a,b,p in [('.4','.7','.6'),('1.1','.3','1.2'),('3','4','.6'),('16.2','16','.8'),('64.4','64','1.2')]:
        a,b,p=map(mp.mpf,[a,b,p]);x=raw_rate(a,b,p);y=tilted(a,b,p);rel=abs(x-y)/y
        rates.append(clean(dict(a=a,b=b,p=p,raw=x,rate=y,relative_difference=rel,passed=rel<mp.mpf('1e-48'))))
    for M in ['.7','4','16','64','256']:
        for t in ['0','.5','1','2']:
            row=return_block(M,'.6','.2',t,24);fine=return_block(M,'.6','.2',t,32)
            diff=abs(row['value']-fine['value'])
            record=clean(dict(coarse=row,fine=fine,difference=diff,rounding_allowance='1e-48',
                passed=diff<=row['tail']+fine['tail']+mp.mpf('1e-48')))
            rows.append(record)
            with raw.open('a',encoding='utf-8') as out:out.write(json.dumps(record)+'\n');out.flush()
            print(M,t,'vacuum',mp.nstr(fine['vacuum_commutator_norm'],8),'passed',record['passed'],flush=True)
    finite=finite_fock();passed=all(r['passed'] for r in rates+rows+finite['bounds']+finite['fourpoint']+finite['commutators'])
    paths=['code/ds3_fock_vacuum.py','code/check_ds3_fock_vacuum.py','code/ds3_two_matter.py','code/ds3_group_pairing.py']
    data=dict(status='passed' if passed else 'failed',rates=rates,vacuum=rows,finite_fock=finite,
        working_dps=60,block_bound=clean(block_bound('.2')),
        dependencies=dict(python=platform.python_version(),mpmath=mp.__version__,numpy=np.__version__),
        scope='Vacuum return and finite Fock-path diagnostics. No fresh crossed U evaluations or all-spin numerical certification. Quadratures are not interval certified.',
        hashes={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in paths})
    dest.write_text(json.dumps(data,indent=2),encoding='utf-8');print(data['status'],flush=True)
    if not passed:raise SystemExit(1)

if __name__=='__main__':main()
