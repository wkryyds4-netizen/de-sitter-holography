"""One-shot complete-plane certificate for the first right F overlap."""
from pathlib import Path
import sys,json,time,hashlib,platform
root=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(root/'tmp/arb-runtime'))
import numpy as np
import flint
from flint import arb,acb,ctx
from chirp_residual_certificate import norm_squared,interval
from tensor_overlap import GaussianOverlap

def main():
    base=root/'results/stage4/residual-certificate'
    output=base/'tensor-gaussian.json';log=base/'tensor-gaussian-cells.jsonl'
    if output.exists() or log.exists():raise RuntimeError('Do not overwrite or restart')
    ctx.prec=256;L=5;terms=3;order=24
    source=root/'results/stage4/flip-domain/coefficients-06.npz'
    with np.load(source,allow_pickle=False) as data:v,w=data['right'][0:2]
    vacuum=np.zeros_like(v);vacuum[0,0]=1;assert np.array_equal(v,vacuum)
    nw=norm_squared(w);normw=nw.sqrt();start=time.perf_counter()
    # Whole complement of the square, with |F|=1.
    tail=(2*arb(L).erfc()*nw).sqrt()
    r=(arb(L)*arb.pi().sqrt()-(2*terms+1)*arb.pi()).exp()
    product=4*r/((1-(-2*arb.pi()).exp())*(1-r))*normw
    evaluator=GaussianOverlap(w,terms,order)
    progress=lambda leaves,e:print(json.dumps(dict(event='partition',leaves=leaves,error=e.str(8))),flush=True)
    cells,error=evaluator.partition(L,arb('0.00002'),progress)
    print(json.dumps(dict(event='partition-complete',leaves=len(cells),
        direct=sum(c['kind']=='direct' for c in cells),error=error.str(15),tail=tail.str(15),product=product.str(15))),flush=True)
    total=acb(0)
    with log.open('x') as f:
        for i,cell in enumerate(cells):
            value=evaluator.quadrature(cell['rect']) if cell['kind']=='gauss' else acb(0)
            total+=value
            row=dict(index=i,rect_dyadic=[list(float(x).as_integer_ratio()) for x in cell['rect']],
                     depth=cell['depth'],kind=cell['kind'],analytic_or_direct_error=interval(cell['error']),
                     value_real=interval(value.real),value_imag=interval(value.imag))
            f.write(json.dumps(row)+'\n');f.flush()
            if i%50==0:print(json.dumps(dict(event='integrate',cell=i,total=len(cells),seconds=time.perf_counter()-start)),flush=True)
    omission=tail+product+error
    real_lower=total.real-omission;deficit=nw-real_lower
    payload=dict(status='full first F overlap enclosed',target_proved=bool(deficit<arb('0.0001')),
        target_eta='0.0001',index=6,n=128,branch='right',step=1,L=L,terms=terms,gauss_order=order,
        precision_bits=ctx.prec,leaf_count=len(cells),direct_cells=sum(c['kind']=='direct' for c in cells),
        input_file=source.relative_to(root).as_posix(),input_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        output_norm_squared=interval(nw),gaussian_tail=interval(tail),product_tail=interval(product),
        cubature_error=interval(error),sum_real=interval(total.real),sum_imag=interval(total.imag),
        real_lower_bound=interval(real_lower),deficit_upper_bound=interval(deficit),seconds=time.perf_counter()-start,
        versions=dict(python=platform.python_version(),python_flint=flint.__version__,numpy=np.__version__),
        hashes={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in
                ['code/tensor_overlap.py','code/run_tensor_gaussian.py','code/chirp_residual_certificate.py',
                 'code/residual_overlap_rectangle.py','results/stage4/residual-certificate/tensor-gaussian-cells.jsonl']},
        scope='Only right step 1 certified. Other F overlaps remain unproved; no complete word counterexample.')
    output.write_text(json.dumps(payload,indent=2))
    print(json.dumps(dict(status=payload['status'],target_proved=payload['target_proved'],
        real=total.real.str(15),deficit=deficit.str(15),seconds=payload['seconds'])),flush=True)

if __name__=='__main__':main()
