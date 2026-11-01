"""One-shot remaining F overlap; choose branch and step explicitly."""
from pathlib import Path
import sys,json,time,hashlib,platform,argparse
root=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(root/'tmp/arb-runtime'))
import numpy as np
import flint
from flint import arb,acb,ctx
from tensor_overlap import GaussianOverlap
from tensor_pair_overlap import PairOverlap,fourier
from chirp_residual_certificate import norm_squared,interval

def main():
    ap=argparse.ArgumentParser();ap.add_argument('branch',choices=['left','right']);ap.add_argument('step',type=int)
    args=ap.parse_args();branch,step=args.branch,args.step
    assert (branch,step) in [('left',1),('left',2),('right',3),('right',5)]
    base=root/'results/stage4/residual-certificate';name=f'tensor-{branch}-{step}'
    output=base/(name+'.json');log=base/(name+'-cells.jsonl');inp=base/(name+'-input.json')
    if any(p.exists() for p in [output,log,inp]):raise RuntimeError('Do not overwrite or restart')
    ctx.prec=256;source=root/'results/stage4/flip-domain/coefficients-06.npz'
    with np.load(source,allow_pickle=False) as data:v,w=data[branch][step-1:step+1]
    transform=(branch,step)!=('left',2)
    if transform:v,w=fourier(v),fourier(w)
    nv,nw=norm_squared(v),norm_squared(w)
    vacuum=np.zeros_like(v);vacuum[0,0]=1;gaussian=np.array_equal(v,vacuum)
    L,terms,order=(5,3,24) if gaussian else (16,8,24)
    ev=GaussianOverlap(w,terms,order) if gaussian else PairOverlap(v,w,terms,order)
    if gaussian:
        tv=2*arb(L).erfc();tw=nw;tail=(tv*tw).sqrt()
    else:tv,tw,tail=ev.tails(L)
    r=(arb(L)*arb.pi().sqrt()-(2*terms+1)*arb.pi()).exp()
    product=4*r/((1-(-2*arb.pi()).exp())*(1-r))*(nv*nw).sqrt()
    inputs=['code/run_tensor_pair.py','code/tensor_pair_overlap.py','code/tensor_overlap.py',
            'code/chirp_residual_certificate.py','code/residual_overlap_rectangle.py',source.relative_to(root).as_posix()]
    meta=dict(branch=branch,step=step,index=6,n=128,fourier_applied=transform,
              gaussian_input=gaussian,L=L,terms=terms,gauss_order=order,precision_bits=ctx.prec,
              input_file=source.relative_to(root).as_posix(),input_norm_squared=interval(nv),output_norm_squared=interval(nw),
              input_tail_mass_bound=interval(tv),output_tail_mass_bound=interval(tw),space_tail=interval(tail),product_tail=interval(product),
              versions=dict(python=platform.python_version(),numpy=np.__version__,python_flint=flint.__version__),
              hashes={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in inputs})
    with inp.open('x') as f:json.dump(meta,f,indent=2)
    start=time.perf_counter()
    cells,error=ev.partition(L,arb('2e-5'))
    print(json.dumps(dict(event='partition',branch=branch,step=step,leaves=len(cells),direct=sum(c['kind']=='direct' for c in cells),tail=tail.str(12),error=error.str(12))),flush=True)
    total=acb(0)
    with log.open('x') as f:
        for i,c in enumerate(cells):
            q=ev.quadrature(c['rect']) if c['kind']=='gauss' else acb(0);total+=q
            row=dict(index=i,rect_dyadic=[list(float(x).as_integer_ratio()) for x in c['rect']],depth=c['depth'],kind=c['kind'],analytic_or_direct_error=interval(c['error']),value_real=interval(q.real),value_imag=interval(q.imag))
            f.write(json.dumps(row)+'\n');f.flush()
            if i%500==0:print(json.dumps(dict(event='integrate',cell=i,total=len(cells),seconds=time.perf_counter()-start)),flush=True)
    deficit=nw-total.real+error+tail+product
    meta.update(status='full-plane F overlap enclosed',target_proved=bool(deficit<arb('1e-4')),target_eta='0.0001',
                cubature_error=interval(error),sum_real=interval(total.real),sum_imag=interval(total.imag),deficit_upper_bound=interval(deficit),
                leaf_count=len(cells),direct_cells=sum(c['kind']=='direct' for c in cells),seconds=time.perf_counter()-start)
    meta['hashes'][log.relative_to(root).as_posix()]=hashlib.sha256(log.read_bytes()).hexdigest()
    with output.open('x') as f:json.dump(meta,f,indent=2)
    print(json.dumps(dict(event='complete',branch=branch,step=step,target_proved=meta['target_proved'],deficit=deficit.str(15),seconds=meta['seconds'])),flush=True)

if __name__=='__main__':main()
