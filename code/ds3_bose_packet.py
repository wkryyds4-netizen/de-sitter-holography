"""Bose 0->2 sector: labelled chains followed by sqrt(2) P_+.

The exchange is the physical matter flip, not a unitary scalar subblock.
Retain raw complex kernels and evaluate the order difference stably.
"""
import argparse
import hashlib
import json
import platform
import time
from pathlib import Path
import numpy as np
import scipy
from ds3_racah_fast import physical
from ds3_crossed_packet import mesh, initial_mesh, g, balanced, packet_weights
from ds3_packet_tail import bounds

ROOT = Path(__file__).resolve().parents[1]


def phase_pair(a, b, c, t):
    ha, hb, hc = (np.sqrt(1+v*v) for v in (a, b, c))
    alpha = np.exp(-1j*t*(ha-hb))
    beta = np.exp(-1j*t*(hb-hc))
    # Exact zero at t=0 and stable for small t.
    delta = 2j*np.exp(-.5j*t*(ha-hc))*np.sin(t*(hb-(ha+hc)/2))
    return alpha, beta, delta


def direct(order, nc, L, times, p=.6, s=.2):
    x,w=mesh(order,s,L);cs,cw=initial_mesh(nc)
    a,b=x[:,None],x[None,:]
    N=0.;Y=np.zeros(len(times),complex);D=np.zeros(len(times))
    moment=0.
    for c,wc in zip(cs,cw):
        wt=wc*w[:,None]*w[None,:]*g(a)/g(c)*balanced(a,b,p)*balanced(b,c,p)
        N+=wt.sum()
        energy=2*np.sqrt(1+b*b)-np.sqrt(1+a*a)-np.sqrt(1+c*c)
        moment+=np.sum(wt*energy**2)
        for it,t in enumerate(times):
            alpha,beta,delta=phase_pair(a,b,c,t)
            Y[it]+=np.sum(wt*np.conj(alpha)*beta)
            D[it]+=np.sum(wt*abs(delta)**2)
    return dict(N=float(N),Y=Y,D=D,moment=float(moment))


def exchange(x,w,cs,cw,U,times,p=.6):
    a,b,d=np.meshgrid(x,x,x,indexing='ij')
    a,b,d=[v.ravel() for v in (a,b,d)]
    w3=(w[:,None,None]*w[None,:,None]*w[None,None,:]).ravel()
    E1=np.zeros(len(times),complex);E2=E1.copy();Z=E1.copy();DX=E1.copy()
    moment=0j;absolute=0.
    for ic,(c,wc) in enumerate(zip(cs,cw)):
        wt=wc*w3*packet_weights(a,b,d,c,p,p)*U[ic]
        absolute+=np.sum(abs(wt))
        eb=2*np.sqrt(1+b*b)-np.sqrt(1+a*a)-np.sqrt(1+c*c)
        ed=2*np.sqrt(1+d*d)-np.sqrt(1+a*a)-np.sqrt(1+c*c)
        moment+=np.sum(wt*eb*ed)
        for it,t in enumerate(times):
            ab,bb,db=phase_pair(a,b,c,t);ad,bd,dd=phase_pair(a,d,c,t)
            E1[it]+=np.sum(wt*np.conj(ab)*ad)
            E2[it]+=np.sum(wt*np.conj(bb)*bd)
            Z[it]+=np.sum(wt*np.conj(ab)*bd)
            DX[it]+=np.sum(wt*np.conj(db)*dd)
    return dict(E1=E1,E2=E2,Z=Z,DX=DX,moment=moment,absolute=float(absolute))


def pairs(v):
    return [[float(z.real),float(z.imag)] for z in np.atleast_1d(v)]


def run():
    parser=argparse.ArgumentParser()
    parser.add_argument('--order',type=int,default=4)
    parser.add_argument('--nc',type=int,default=8)
    parser.add_argument('--L',type=float,default=8)
    parser.add_argument('--series',type=int,default=128)
    parser.add_argument('--tail-order',type=int,default=12)
    parser.add_argument('--label',required=True)
    args=parser.parse_args()
    out=ROOT/'results/stage4/ds3-bose';out.mkdir(parents=True,exist_ok=True)
    dest=out/(args.label+'.json');raw=out/(args.label+'.npz')
    if dest.exists() or raw.exists():raise RuntimeError('Use a fresh run label')
    x,w=mesh(args.order,L=args.L);cs,cw=initial_mesh(args.nc)
    a,b,d=np.meshgrid(x,x,x,indexing='ij');a,b,d=[v.ravel() for v in (a,b,d)]
    U=np.empty((len(cs),len(a)),complex);cond=np.empty(U.shape)
    start=time.time()
    data=dict(status='running',parameters=vars(args),model=dict(s=.2,p=.6,G=0,radius=1),
              times=[0.,.001,.01,.5,1.,2.],
              state='rho|f|^2=Beta(5,5) on [.4,.8], real nonnegative f; half-period creators',
              versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__))
    for ic,c in enumerate(cs):
        for lo in range(0,len(a),256):
            sl=slice(lo,min(lo+256,len(a)))
            U[ic,sl],cond[ic,sl]=physical(b[sl],d[sl],[a[sl],.6,c,.6],args.series,args.tail_order)
        if not np.all(np.isfinite(U[ic])):raise RuntimeError('Nonfinite kernel')
        data.update(completed_initial_nodes=ic+1,total_initial_nodes=len(cs),seconds=time.time()-start)
        dest.write_text(json.dumps(data,indent=2),encoding='utf-8')
        print(args.label,ic+1,len(cs),round(time.time()-start,2),flush=True)
    np.savez_compressed(raw,x=x,w=w,c=cs,cw=cw,U=U,condition=cond)
    ts=data['times'];di=direct(args.order,args.nc,args.L,ts);ex=exchange(x,w,cs,cw,U,ts)
    hi=direct(32,16,20,ts)
    bnd=bounds(args.L,p='.6',q='.6')
    eps=float(bnd['tail_norm_bound']);nmax=float(bnd['norm_squared_bound'])
    rows=[]
    for j,t in enumerate(ts):
        n1=di['N']+ex['E1'][j];n2=di['N']+ex['E2'][j]
        ov=di['Y'][j]+ex['Z'][j]
        stable=di['D'][j]+ex['DX'][j]
        rows.append(dict(t=t,N_label=di['N'],E1=pairs(ex['E1'][j])[0],E2=pairs(ex['E2'][j])[0],
                         Y=pairs(di['Y'][j])[0],Z=pairs(ex['Z'][j])[0],
                         N1=pairs(n1)[0],N2=pairs(n2)[0],overlap=pairs(ov)[0],
                         D_direct=float(di['D'][j]),D_exchange=pairs(ex['DX'][j])[0],
                         D=pairs(stable)[0],expanded_D=float((n1+n2-2*ov.real).real),
                         cauchy_ratio=float(abs(ov)/np.sqrt(n1.real*n2.real))))
    paths=['code/ds3_bose_packet.py','code/ds3_racah_fast.py','code/ds3_crossed_packet.py','code/ds3_packet_tail.py']
    data.update(status='completed',seconds=time.time()-start,rows=rows,points=int(U.size),
                raw_grid=raw.relative_to(ROOT).as_posix(),raw_grid_sha256=hashlib.sha256(raw.read_bytes()).hexdigest(),
                hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
                direct_high=dict(N=hi['N'],Y=pairs(hi['Y']),D=hi['D'].tolist(),moment=hi['moment'],order=32,nc=16,L=20),
                moment_direct=di['moment'],moment_exchange=pairs(ex['moment'])[0],
                moment_total=pairs(di['moment']+ex['moment'])[0],absolute_sum=ex['absolute'],
                max_condition=float(cond.max()),max_imaginary_kernel=float(abs(U.imag).max()),
                mass_tail=dict(epsilon=eps,Nmax=nmax,D_bound=4*eps+8*np.sqrt(eps*nmax),
                               bose_inner_product_bound=2*eps+2*np.sqrt(eps*nmax)),
                error_scope='Mass tails analytic; quadrature, special-function and floating-point errors empirical.')
    # General Bose inner product has direct and exchange tails: epsilon + 2 sqrt(epsilon Nmax).
    data['mass_tail']['bose_inner_product_bound']=eps+2*np.sqrt(eps*nmax)
    dest.write_text(json.dumps(data,indent=2),encoding='utf-8')
    print(json.dumps(dict(rows=rows,moment=data['moment_total'],mass_tail=data['mass_tail']),indent=2),flush=True)


if __name__=='__main__':run()
