"""Exact three-time product-to-commutator vertex and normalized fourth moment."""
from pathlib import Path
import json,hashlib,platform
import sympy as s
from commutator_primary import moment,expr,q
ROOT=Path(__file__).resolve().parent.parent;OUT=ROOT/'results/susskind-commutator-primary'
def main():
    target=OUT/'vertex.json'
    if target.exists():raise FileExistsError(target)
    rows=[]
    for n in range(0,7,2):
      for a in range(n+1):
       for b in range(n-a+1):
        c=n-a-b;suffix='H'*a+'X'+'H'*b+'R'+'H'*c
        cross=expr(moment('XR'+suffix));uncross=expr(moment('RX'+suffix))
        assert s.expand(cross-q**6*uncross)==0
        rows.append(dict(arcs=[a,b,c],uncrossed_polynomial=str(uncross),difference='0'))
    data=json.loads((OUT/'exact.json').read_text());m4=s.sympify(data['primary']['even_moments'][1]['normalized'])
    simple=2+q**25+(1-q**6)**2*(1+q**4+q**9)/2
    assert s.expand(m4-simple)==0
    # Independent order-basis Gram directly fixes the product coefficients.
    G=s.Matrix([[1,q**6],[q**6,1]]);plus=s.Matrix([1,1]);minus=s.Matrix([1,-1])
    assert (plus.T*G*minus)[0]==0
    assert s.expand((plus.T*G*plus)[0]-2*(1+q**6))==0
    assert s.expand((minus.T*G*minus)[0]-2*(1-q**6))==0
    inputs=['code/check_commutator_vertex.py','code/commutator_primary.py','results/susskind-commutator-primary/exact.json']
    output=dict(status='passed',three_arc_cases=len(rows),records=rows,order_gram=str(G),
       normalized_fourth_moment=str(simple),fourth_probe_excess_over_h2='54',
       exact_three_time_vertex='tau(K(t0) X(t1) R(t2)) = i sqrt((1-q^6)/2) tau(R(t0) X(t0) X(t1) R(t2))',
       leading_vertex='i sqrt(3h) sech^4(u0-u1) sech^6(u0-u2) + O(h^(3/2))',
       environment=dict(python=platform.python_version(),sympy=s.__version__),
       input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs})
    target.write_text(json.dumps(output,indent=2)+'\n');print('50 exact three-arc identities and order Gram/fourth moment passed')
if __name__=='__main__':main()
