"""Exact-rational coverage/accounting and ball-bound replay of first F step."""
from pathlib import Path
from fractions import Fraction as Q
import sys,json,hashlib
root=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(root/'tmp/arb-runtime'))
import numpy as np
from flint import arb,ctx
from tensor_overlap import GaussianOverlap

def ends(item):return tuple(Q(int(a),int(b)) for a,b in [item['lower'],item['upper']])
def overlap(a,b):return a[0]<=b[1] and b[0]<=a[1]
def ball(x):return (Q(str(x.lower().fmpq())),Q(str(x.upper().fmpq())))
def main():
    base=root/'results/stage4/residual-certificate'
    p=json.loads((base/'tensor-gaussian.json').read_text())
    rows=[json.loads(s) for s in (base/'tensor-gaussian-cells.jsonl').read_text().splitlines()]
    checks=[]
    def check(name,ok):checks.append(dict(check=name,passed=bool(ok)))
    for path,h in p['hashes'].items():check('hash '+path,hashlib.sha256((root/path).read_bytes()).hexdigest()==h)
    check('input hash',hashlib.sha256((root/p['input_file']).read_bytes()).hexdigest()==p['input_sha256'])
    paths=set();area=Q(0)
    for row in rows:
        a,b,c,d=[Q(*x) for x in row['rect_dyadic']]
        depth=row['depth'];size=Q(10,2**depth)
        ix=(a+5)/size;iy=(c+5)/size
        check('dyadic square '+str(row['index']),b-a==d-c==size and ix.denominator==iy.denominator==1 and 0<=ix<2**depth and 0<=iy<2**depth)
        path=tuple(((int(ix)>>k)&1)+2*((int(iy)>>k)&1) for k in reversed(range(depth)))
        check('unique leaf '+str(row['index']),path not in paths);paths.add(path);area+=(b-a)*(d-c)
    check('no overlapping ancestors',not any(path[:k] in paths for path in paths for k in range(len(path))))
    check('exact complete square area',area==100)
    check('leaf count',len(rows)==p['leaf_count'])
    e=sum((ends(r['analytic_or_direct_error'])[1] for r in rows),Q(0))
    err=ends(p['cubature_error']);check('stored error encloses sum',err[0]<=e<=err[1])
    sums={part:(sum((ends(r['value_'+part])[0] for r in rows),Q(0)),sum((ends(r['value_'+part])[1] for r in rows),Q(0))) for part in ['real','imag']}
    for part in sums:check('sum '+part,overlap(sums[part],ends(p['sum_'+part])))
    omitted=e+ends(p['gaussian_tail'])[1]+ends(p['product_tail'])[1]
    deficit=ends(p['output_norm_squared'])[1]-sums['real'][0]+omitted
    check('target from exact rational worst endpoints',deficit<Q(1,10000))
    ctx.prec=320
    with np.load(root/p['input_file'],allow_pickle=False) as data:w=data['right'][1]
    ev=GaussianOverlap(w,p['terms'],p['gauss_order'])
    for row in rows:
        rect=[float(Q(*x)) for x in row['rect_dyadic']]
        new=ev.bound(tuple(rect),row['depth'])
        check('bound replay '+str(row['index']),new['kind']==row['kind'] and overlap(ball(new['error']),ends(row['analytic_or_direct_error'])))
    for index in [0,1,100,500,924]:
        row=rows[index]
        if row['kind']!='gauss':continue
        value=ev.quadrature(tuple(float(Q(*x)) for x in row['rect_dyadic']))
        for part in ['real','imag']:check('higher precision quadrature '+str(index)+' '+part,overlap(ball(getattr(value,part)),ends(row['value_'+part])))
    out=dict(passed=all(c['passed'] for c in checks),count=len(checks),checks=checks,
             exact_deficit_upper=[str(deficit.numerator),str(deficit.denominator)],
             scope='First right F overlap only; shared library and mathematical bounds are not independently formally verified.',
             hashes={path:hashlib.sha256((root/path).read_bytes()).hexdigest() for path in
                     ['code/audit_tensor_gaussian.py','results/stage4/residual-certificate/tensor-gaussian.json']})
    with (base/'tensor-gaussian-audit.json').open('x') as f:json.dump(out,f,indent=2)
    print(json.dumps(dict(passed=out['passed'],count=out['count'],failed=[c for c in checks if not c['passed']])))

if __name__=='__main__':main()
