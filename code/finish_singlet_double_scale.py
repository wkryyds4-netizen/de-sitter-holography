"""Finish the local stage after the isolated sparse-overlap audit repair."""
from pathlib import Path
from fractions import Fraction
import hashlib,json,platform
import numpy as np
import sympy as sp
from singlet_coefficient import overlap_table,local_polynomial,pairings,moment_coefficients

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-singlet-double-scale'


def save(name,obj):
    with (OUT/name).open('x',encoding='utf-8') as f:json.dump(obj,f,indent=2);f.write('\n')


def main():
    for name in ['overlap-graphs.json','validation.json']:
        if (OUT/name).exists():raise FileExistsError(name)
    graph=overlap_table();grouped={c:[0,0] for c in range(4)}
    for r in graph:
        grouped[r['crossings']][0]+=r['sign']*r['delta']
        grouped[r['crossings']][1]+=r['sign']*r['alpha']*r['beta']
    assert grouped=={0:[8,0],1:[0,32],2:[-8,-32],3:[0,0]}
    # All Hamiltonian sixth-moment diagrams must recover the known q chord factor.
    for pair in pairings(tuple(range(6))):
        terms=local_polynomial(6,pair,((-2,-1),(-4,-3)))
        single={b.index(1):v[0] for b,v in terms if sum(b)==1}
        kappa=-15
        for bits,v in terms:
            if sum(bits)==2:
                ij=[i for i,b in enumerate(bits) if b]
                kappa+=v[0]*single[ij[0]]*single[ij[1]]
        cr=sum(a<c<b<d or c<a<d<b for i,(a,b) in enumerate(pair) for c,d in pair[i+1:])
        assert kappa-3==-4*cr
    saved=json.loads((OUT/'coefficients-p2.json').read_text())['coefficients']
    assert moment_coefficients(2)==moment_coefficients(2,'all15')==saved
    tensor=json.loads((OUT/'independent-tensor.json').read_text())
    assert Fraction(tensor['total']['mu2'])==2700
    assert tensor['integer_contraction_max']==0
    words=json.loads((OUT/'word-checks.json').read_text())
    assert words['scalar_comparisons']==2880 and words['status']=='passed'
    oldcenter=json.loads((OUT/'centered-gram-initial.json').read_text())
    newcenter=json.loads((OUT/'centered-gram.json').read_text())
    assert oldcenter['rows']==newcenter['rows']
    save('overlap-graphs.json',dict(rows=graph,lambda_and_lambda_squared_coefficients=grouped,
                                   unmarked_sixth_moment_diagrams_checked=15))
    files=['code/singlet_coefficient.py','code/check_singlet_coefficient.py','code/run_singlet_double_scale.py',
       'code/run_singlet_centered_moment.py','code/finish_singlet_double_scale.py',
       'code/singlet_wick_moments.py','code/complex_syk_dictionary.py',
       'results/susskind-singlet-moments/moments.json']
    save('validation.json',dict(status='passed after recorded audit repair',
       environment=dict(python=platform.python_version(),numpy=np.__version__,sympy=sp.__version__),
       input_sha256={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in files},
       exact_checks=dict(marked_words=720,scalar_comparisons=2880,
          nonzero_mixed_words=words['nonzero_mixed_words'],negative_mixed_words=words['negative_mixed_words'],
          finite_coefficients_p=[1,2,3,4,5],all15_cross_equal_p=[1,2,3],
          independent_full_tensor_N5_p2='2700',overlap_graphs=24,known_chord_weights=15),
       repair='Missing -1 contributions for zero minors in exponent audit; initial code and record retained. Coefficient extraction and all computed moments were unchanged.',
       scope='Finite coefficients exact; double-scaled limit from the written overlap proof. Finite-p ratios have no claimed asymptotic error bar.'))
    print(json.dumps(dict(status='passed',word_cases=720,independent_moment='2700',graph_groups=grouped)))


if __name__=='__main__':main()
