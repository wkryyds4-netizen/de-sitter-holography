"""Checkpointed physical packet quadrature, retaining complete complex kernels."""
import argparse
import hashlib
import json
import platform
import time
from pathlib import Path
import numpy as np
import scipy
from ds3_racah_fast import physical
from ds3_crossed_packet import mesh,initial_mesh,packet_weights,positive_norms

ROOT=Path(__file__).resolve().parents[1]


def run():
    parser=argparse.ArgumentParser()
    parser.add_argument('--order',type=int,default=4)
    parser.add_argument('--nc',type=int,default=8)
    parser.add_argument('--L',type=float,default=8)
    parser.add_argument('--series',type=int,default=128)
    parser.add_argument('--tail-order',type=int,default=12)
    parser.add_argument('--label',required=True)
    args=parser.parse_args()
    out=ROOT/'results/stage4/ds3-crossed-packet';out.mkdir(parents=True,exist_ok=True)
    dest=out/(args.label+'.json');gridpath=out/(args.label+'.npz')
    if dest.exists() or gridpath.exists():raise RuntimeError('Preserve prior run; choose a new label')
    x,w=mesh(args.order,L=args.L);cs,cw=initial_mesh(args.nc)
    a,b,d=np.meshgrid(x,x,x,indexing='ij');a,b,d=[v.ravel() for v in [a,b,d]]
    w3=(w[:,None,None]*w[None,:,None]*w[None,None,:]).ravel()
    U=np.empty((len(cs),len(a)),complex);cond=np.empty(U.shape,float)
    times=np.array([0.,.5,1.,2.]);values=np.zeros(len(times),complex)
    absolute=0.;proxy=0.;max_imag=0.;worst=[];start=time.time()
    paths=['code/ds3_racah_fast.py','code/ds3_crossed_packet.py','code/run_ds3_crossed_packet.py']
    data=dict(parameters=vars(args),state='rho(c)|f(c)|^2=Beta(5,5) on [.4,.8]; f real nonnegative',
              model=dict(radius=1,G=0,s=.2,p=.6,q=.9,half_period=True),
              versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
              hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},status='running')
    for ic,(c,wc) in enumerate(zip(cs,cw)):
        for lo in range(0,len(a),256):
            sl=slice(lo,min(lo+256,len(a)))
            z,k=physical(b[sl],d[sl],[a[sl],.9,c,.6],args.series,args.tail_order)
            if not np.all(np.isfinite(z)):raise RuntimeError('Nonfinite kernel requires a separate precision repair')
            U[ic,sl]=z;cond[ic,sl]=k
            wt=wc*w3[sl]*packet_weights(a[sl],b[sl],d[sl],c)
            phase_energy=np.sqrt(1+a[sl]**2)+np.sqrt(1+c*c)-np.sqrt(1+b[sl]**2)-np.sqrt(1+d[sl]**2)
            values+=np.sum((wt*z)[:,None]*np.exp(1j*phase_energy[:,None]*times),axis=0)
            absolute+=float(np.sum(wt*abs(z)))
            sensitivity=wt*abs(z)*k*np.finfo(float).eps
            proxy+=float(np.sum(sensitivity));max_imag=max(max_imag,float(max(abs(z.imag))))
            j=int(np.argmax(sensitivity));worst.append(dict(a=float(a[sl][j]),b=float(b[sl][j]),d=float(d[sl][j]),c=float(c),
                u=[float(z[j].real),float(z[j].imag)],condition=float(k[j]),weighted_roundoff_proxy=float(sensitivity[j])))
        data.update(completed_initial_nodes=ic+1,total_initial_nodes=len(cs),seconds=time.time()-start)
        dest.write_text(json.dumps(data,indent=2),encoding='utf-8')
        print(args.label,ic+1,len(cs),'seconds',round(time.time()-start,2),flush=True)
    np.savez_compressed(gridpath,x=x,w=w,c=cs,cw=cw,U=U,condition=cond)
    norms=positive_norms(order=32,nc=16,L=20)
    data.update(status='completed',seconds=time.time()-start,points=int(U.size),
        times=times.tolist(),values=[[float(z.real),float(z.imag)] for z in values],
        absolute_sum=absolute,cancellation=[float(absolute/abs(v)) for v in values],
        roundoff_proxy=proxy,roundoff_proxy_is_rigorous=False,max_condition=float(np.max(cond)),
        max_kernel_imag=max_imag,worst_weighted_nodes=sorted(worst,key=lambda r:r['weighted_roundoff_proxy'],reverse=True)[:12],
        positive_norms=norms,raw_grid=gridpath.relative_to(ROOT).as_posix(),raw_grid_sha256=hashlib.sha256(gridpath.read_bytes()).hexdigest())
    dest.write_text(json.dumps(data,indent=2),encoding='utf-8')
    print(json.dumps({k:data[k] for k in ['values','roundoff_proxy','max_condition','positive_norms','seconds']},indent=2),flush=True)


if __name__=='__main__':run()
