"""Exact primitive and one-selected density identities in the sealed transfer."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from datetime import datetime, timezone
import hashlib,json,platform
import sympy as sp
from slater_covariant_transfer import step,inner,moment,add,scale

ROOT=Path(__file__).resolve().parent.parent
DIR=ROOT/'results/susskind-slater-neutral-slope';OUT=DIR/'amplification-checks.json'
if OUT.exists():raise FileExistsError(OUT)
rows=[];density_count=0;primitive_count=0;tower_count=0
for gamma,d in ((F(1,2),F(1,2)),(F(2,3),F(1,2))):
    c=-d*gamma; r=1-gamma**2-c**2
    zeta={'P':F(1),'H':-gamma,'U':-c}
    assert not step(zeta,'a',gamma,d) and not step(zeta,'v',gamma,d)
    assert inner(zeta,zeta,gamma,d)==r>0
    for length in range(9):
        for letters in product('HU',repeat=length):
            word=''.join(letters)
            rhs=r*moment(word,gamma,d)
            for a,ca in (('H',gamma),('U',c)):
                for b,cb in (('H',gamma),('U',c)):
                    rhs+=ca*cb*moment(a+word+b,gamma,d)
            assert moment('-'+word+'+',gamma,d)==rhs,word
            density_count+=1
    primitive_rows=[]
    for level in (1,2):
        words=[''.join(x) for x in product('HUP',repeat=level)]
        lower=[''.join(x) for x in product('HUP',repeat=level-1)]
        mat=sp.Matrix([[sp.Rational(step({w:F(1)},op,gamma,d).get(v,0))
                        for w in words] for op in ('a','v') for v in lower])
        kernels=mat.nullspace(); prims=[]
        for col in kernels:
            vec={w:F(x) for w,x in zip(words,col) if x}
            assert not step(vec,'a',gamma,d) and not step(vec,'v',gamma,d)
            # Gram-Schmidt only in the actual transfer metric.
            for prev in prims:
                vec=add(vec,scale(-inner(prev,vec,gamma,d)/inner(prev,prev,gamma,d),prev))
            assert inner(vec,vec,gamma,d)>0
            prims.append(vec);primitive_count+=1
        for vec in prims:
            norm=inner(vec,vec,gamma,d)
            for a in ('H','U'):
                for b in ('H','U'):
                    aa=step(vec,'A' if a=='H' else 'V',gamma,d)
                    bb=step(vec,'A' if b=='H' else 'V',gamma,d)
                    assert inner(aa,bb,gamma,d)==norm*(a==b)
                    tower_count+=1
        primitive_rows.append(dict(level=level,kernel_dimension=len(prims),
                                   primitive_squared_norms=[str(inner(v,v,gamma,d)) for v in prims]))
    rows.append(dict(gamma=str(gamma),sqrt_lambda=str(d),new_joint_primitive_norm=str(r),
                     primitive_checks=primitive_rows))
DIR.mkdir(parents=True,exist_ok=True)
files=[Path(__file__),ROOT/'code/slater_covariant_transfer.py',ROOT/'results/susskind-slater-covariant/local-evidence.json']
result=dict(status='passed',created_at_utc=datetime.now(timezone.utc).isoformat(),rows=rows,
    counts=dict(density_identities=density_count,primitive_vectors=primitive_count,tower_inner_products=tower_count),
    parameters=dict(rational_pairs=[['1/2','1/2'],['2/3','1/2']],max_density_word_degree=8,primitive_levels=[1,2],physical_relation_imposed=False),
    environment=dict(python=platform.python_version(),sympy=sp.__version__,arithmetic='fractions.Fraction and SymPy Rational'),
    input_sha256={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
    scope='Finite exact checks of the graded q-Fock amplification proof and selected density formula; not an extrapolation of finite kernel ranks.')
with OUT.open('x',encoding='utf-8') as stream:json.dump(result,stream,indent=2);stream.write('\n')
print(json.dumps(dict(status='passed',counts=result['counts'],sha256=hashlib.sha256(OUT.read_bytes()).hexdigest())))
