"""Reconstruct actual quadratic operators from the saved rank filters."""
from fractions import Fraction as F
import hashlib,json,platform
import numpy as np
import sympy as s
from memory_quadratic_contraction import Contraction,anti
from memory_rank_blocks import ROOT,OUT


def exact_matrix(x):
    return s.Matrix(x).applyfunc(lambda z:s.nsimplify(z,rational=True))


def main():
    path=OUT/'rank-reconstruction.json'
    if path.exists():raise FileExistsError(path)
    data=json.loads((OUT/'rank-blocks.json').read_text());o=Contraction(data['N'],data['p'],data['m'])
    pin=[(np.array(r['numerator'],dtype=np.int64),r['denominator']) for r in data['input_projectors']]
    pout=[(s.Matrix(r['numerator']),r['denominator']) for r in data['output_projectors']]
    table=[]
    for row in data['rows']:
        G=s.Matrix([[s.Rational(x) for x in rr] for rr in row['gram']])
        if G[0,0]==0:
            assert G==s.zeros(4);continue
        ratios=[s.cancel(G[0,i]/G[0,0]) for i in range(4)]
        assert G==G[0,0]*s.Matrix(ratios)*s.Matrix(ratios).T
        table.append(dict(input_ranks=row['input_ranks'],output_rank=row['output_rank'],
                          product_weight=str(G[0,0]),multipliers=list(map(str,ratios))))
    reconstruct=[]
    for case in range(2):
        x=np.array([((7+2*case)*i+3+case)%(11+2*case)-5 for i in range(len(o.rs))],dtype=np.int64)
        h=np.einsum('a,aij->ij',x,o.hs);k=np.einsum('a,aij->ij',x,o.ks)
        t=o.tlinear(h,k)+86*h
        actual=[h@h,o.A(h,k),anti(h,t)/2,o.contract(h,k)]
        split=[]
        for num,den in pin:
            xx=num@x;split.append(exact_matrix(np.einsum('a,aij->ij',xx,o.hs))/den)
        assert sum(split,s.zeros(o.d))==exact_matrix(h)
        totals=[s.zeros(o.d) for _ in range(4)]
        for row in table:
            a,b=row['input_ranks'];r=row['output_rank']
            prod=split[a]*split[b]
            if a!=b:prod+=split[b]*split[a]
            pn,pd=pout[r];filtered=(pn*s.Matrix(list(prod))/pd).reshape(o.d,o.d)
            for i,g in enumerate(row['multipliers']):totals[i]+=s.Rational(g)*filtered
        for predicted,true in zip(totals,actual):assert predicted==exact_matrix(true)
        reconstruct.append(dict(case=case,integer_couplings=x.tolist(),operator_identities=4))
    # Full tensor equality follows from every block's exact rank-one Gram;
    # the separate realization reconstruction checks the rank-filter convention.
    names=['code/check_memory_rank_reconstruction.py','code/memory_rank_blocks.py',
           'code/memory_quadratic_contraction.py','results/susskind-memory-rank-blocks/rank-blocks.json']
    result=dict(status='passed',filter_table=table,finite_reconstructions=reconstruct,
        scope='All eleven nonzero blocks of four tested maps are collinear. This does not say the representation tensor products are multiplicity-free, nor that general N,p,m have this property.',
        environment=dict(python=platform.python_version(),numpy=np.__version__,sympy=s.__version__),
        input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names})
    path.write_text(json.dumps(result,indent=2)+'\n')
    print('Exact filters:',len(table),'blocks; independent operator identities:',8,flush=True)


if __name__=='__main__':main()
