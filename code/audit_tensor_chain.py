"""Exact coverage/accounting, bound replay, and a posteriori word estimate.

Only the explicitly defined hbar=1 L2(R2) word is tested. This does not
identify a physical observer operator, gluing rule, or full dS dual.
"""
from pathlib import Path
from fractions import Fraction as Q
import sys,json,hashlib
root=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(root/'tmp/arb-runtime'))
import numpy as np
from flint import arb,ctx,fmpq
from tensor_overlap import GaussianOverlap
from tensor_pair_overlap import PairOverlap,fourier
from chirp_residual_certificate import interval

def ends(item):return tuple(Q(int(a),int(b)) for a,b in [item['lower'],item['upper']])
def overlap(a,b):return a[0]<=b[1] and b[0]<=a[1]
def ball(x):return Q(str(x.lower().fmpq())),Q(str(x.upper().fmpq()))
def lift(q):return arb(fmpq(q.numerator,q.denominator))
def exact_norm(c):return sum((Q(float(z.real))**2+Q(float(z.imag))**2 for z in c.flat),Q(0))

def audit(root=root,replay_bounds=True):
    base=root/'results/stage4/residual-certificate'
    count=0;failed=[];records=[];hashes={}
    def check(name,ok):
        nonlocal count
        count+=1
        if not ok:failed.append(name)
    with np.load(root/'results/stage4/flip-domain/coefficients-06.npz',allow_pickle=False) as data:
        left,right=data['left'],data['right']
    specs=[('tensor-gaussian','right',1),('tensor-left-1','left',1),('tensor-left-2','left',2),('tensor-right-3','right',3),('tensor-right-5','right',5)]
    for name,branch,step in specs:
        p=json.loads((base/(name+'.json')).read_text())
        hashes[(base/(name+'.json')).relative_to(root).as_posix()]=hashlib.sha256((base/(name+'.json')).read_bytes()).hexdigest()
        rows=[json.loads(s) for s in (base/(name+'-cells.jsonl')).read_text().splitlines()]
        for path,h in p['hashes'].items():check(name+' hash '+path,hashlib.sha256((root/path).read_bytes()).hexdigest()==h)
        vec=left if branch=='left' else right;v,w=vec[step-1:step+1]
        nv,nw=exact_norm(v),exact_norm(w)
        check(name+' norm w',ends(p['output_norm_squared'])[0]<=nw<=ends(p['output_norm_squared'])[1])
        if name!='tensor-gaussian':check(name+' norm v',ends(p['input_norm_squared'])[0]<=nv<=ends(p['input_norm_squared'])[1])
        transform=(branch,step) in [('left',1),('right',3),('right',5)]
        if transform:
            v,w=fourier(v),fourier(w)
            check(name+' exact Fourier norm',exact_norm(v)==nv and exact_norm(w)==nw)
        paths=set();area=Q(0);L=p['L']
        for row in rows:
            a,b,c,d=[Q(*x) for x in row['rect_dyadic']]
            depth=row['depth'];size=Q(2*L,2**depth);ix=(a+L)/size;iy=(c+L)/size
            check(name+' square '+str(row['index']),b-a==d-c==size and ix.denominator==iy.denominator==1 and 0<=ix<2**depth and 0<=iy<2**depth)
            path=tuple(((int(ix)>>k)&1)+2*((int(iy)>>k)&1) for k in reversed(range(depth)))
            check(name+' unique '+str(row['index']),path not in paths);paths.add(path);area+=(b-a)*(d-c)
        check(name+' no overlapping ancestors',not any(path[:k] in paths for path in paths for k in range(len(path))))
        check(name+' exact coverage',area==4*L*L)
        check(name+' leaf count',len(rows)==p['leaf_count'])
        e=sum((ends(r['analytic_or_direct_error'])[1] for r in rows),Q(0))
        check(name+' error sum',ends(p['cubature_error'])[0]<=e<=ends(p['cubature_error'])[1])
        sums={part:(sum((ends(r['value_'+part])[0] for r in rows),Q(0)),sum((ends(r['value_'+part])[1] for r in rows),Q(0))) for part in ['real','imag']}
        for part in sums:check(name+' sum '+part,overlap(sums[part],ends(p['sum_'+part])))
        tail_key='gaussian_tail' if name=='tensor-gaussian' else 'space_tail'
        omitted=e+ends(p[tail_key])[1]+ends(p['product_tail'])[1]
        real_lower=sums['real'][0]-omitted;deficit=nw-real_lower;e2=nv+nw-2*real_lower
        check(name+' target',deficit<Q(1,10000));check(name+' nonnegative squared residual bound',e2>=0)
        ctx.prec=256;vacuum=np.zeros_like(v);vacuum[0,0]=1;gaussian=np.array_equal(v,vacuum)
        def evaluator():return GaussianOverlap(w,p['terms'],p['gauss_order']) if gaussian else PairOverlap(v,w,p['terms'],p['gauss_order'])
        ev=evaluator()
        if not gaussian:
            tv,tw,tail=ev.tails(L)
            check(name+' tail replay',overlap(ball(tail),ends(p[tail_key])))
        if replay_bounds:
            for row in rows:
                rect=tuple(float(Q(*x)) for x in row['rect_dyadic'])
                new=ev.bound(rect,row['depth'])
                check(name+' bound '+str(row['index']),new['kind']==row['kind'] and overlap(ball(new['error']),ends(row['analytic_or_direct_error'])))
            ctx.prec=320;ev=evaluator()
            for index in [0,1,len(rows)//3,len(rows)//2,len(rows)-1]:
                row=rows[index]
                if row['kind']!='gauss':continue
                q=ev.quadrature(tuple(float(Q(*x)) for x in row['rect_dyadic']))
                for part in ['real','imag']:check(name+' high precision '+str(index)+' '+part,overlap(ball(getattr(q,part)),ends(row['value_'+part])))
        ctx.prec=512;eb=lift(e2).sqrt()
        records.append(dict(branch=branch,step=step,deficit_upper=[str(deficit.numerator),str(deficit.denominator)],residual_upper=interval(eb.upper()),leaves=len(rows)))
        print(json.dumps(dict(event='audited',branch=branch,step=step,errors=len(failed),residual=eb.str(14))),flush=True)
    chirp_path=base/'chirp.json';chirp=json.loads(chirp_path.read_text())
    c=next(x for x in chirp['rows'] if x['index']==6 and x['precision_bits']==1024)
    check('chirp input hash',hashlib.sha256((root/c['input_file']).read_bytes()).hexdigest()==c['input_sha256'])
    hashes[chirp_path.relative_to(root).as_posix()]=hashlib.sha256(chirp_path.read_bytes()).hexdigest()
    separation2=sum(((Q(float(a.real))-Q(float(b.real)))**2+(Q(float(a.imag))-Q(float(b.imag)))**2 for a,b in zip(left[-1].flat,right[-1].flat)),Q(0))
    ctx.prec=512;separation=lift(separation2).sqrt()
    total=sum((lift(ends(r['residual_upper'])[1]) for r in records),arb(0))
    for row in c['chirp_steps']:total+=lift(ends(row['residual'])[1])
    margin=separation-total
    check('positive continuous word bound',margin>0)
    return dict(passed=not failed,checks=count,failed=failed,overlaps=records,separation=interval(separation),total_step_error_upper=interval(total),word_norm_lower_bound=interval(margin.lower()),
                scientific_status='Specified hbar=1 bounded L2(R2) word has a positive norm lower bound' if margin>0 and not failed else 'No certified conclusion',
                scope='Uses existing quadrature and Mehler bounds with ball arithmetic. Not formal verification; not a physical observer, finite-G or full holographic result. Identification with any published theorem requires its separate representation/definition audit.',hashes=hashes)

def main():
    out=root/'results/stage4/residual-certificate/tensor-chain-audit.json'
    if out.exists():raise RuntimeError('Preserve existing audit')
    result=audit();result['hashes']['code/audit_tensor_chain.py']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    with out.open('x') as f:json.dump(result,f,indent=2)
    print(json.dumps(dict(passed=result['passed'],checks=result['checks'],failed=result['failed'],margin=result['word_norm_lower_bound']['display'])))

if __name__=='__main__':main()
