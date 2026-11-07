"""Exact q-Fock / independent colored-chord tests of four-copy regression.

New output only. These identities validate the limiting algebra, not finite-N
strong convergence. All arithmetic in this driver is rational.
"""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
from collections import defaultdict
from itertools import product
import datetime
import hashlib
import json
import platform
import time

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / 'results/susskind-memory-two-component'


def add(a, b, scale=F(1)):
    out = list(a) + [F(0)] * max(0, len(b)-len(a))
    for i, value in enumerate(b):
        out[i] += scale * value
    return out


def mul(a, b):
    out = [F(0)] * (len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return out


def hermites(q, C, degree):
    polys = [[F(1)], [F(0), F(1)]]
    for n in range(1, degree):
        polys.append(add([F(0)]+polys[-1], polys[-2],
                         -C*sum(q**j for j in range(n))))
    return polys[:degree+1]


def basis_coefficients(poly, basis):
    rem = list(poly)
    out = [F(0)] * len(poly)
    for n in range(len(poly)-1, -1, -1):
        out[n] = rem[n]
        rem = add(rem, basis[n], -out[n])
    assert not any(rem)
    return out


def fock_step(state, color, q, C):
    out = defaultdict(F)
    for word, coeff in state.items():
        out[(color,)+word] += coeff
        for j, letter in enumerate(word):
            if letter == color:
                out[word[:j]+word[j+1:]] += coeff*C*q**j
    return {word: coeff for word, coeff in out.items() if coeff}


def word_state(word, q, C):
    state = {(): F(1)}
    for color in reversed(word):
        state = fock_step(state, color, q, C)
    return state


def conditional(word, q, C, basis):
    out = [F(0)]
    for tensor, coeff in word_state(word, q, C).items():
        if all(x == 'h' for x in tensor):
            out = add(out, basis[len(tensor)], coeff)
    return out


@lru_cache(None)
def colored_pairings(word):
    """Enumerate by matching colors; independent of the Fock recursion."""
    def recurse(points):
        if not points:
            yield ()
            return
        first = points[0]
        for j, second in enumerate(points[1:], 1):
            if word[first] != word[second]:
                continue
            for rest in recurse(points[1:j]+points[j+1:]):
                yield ((first, second),)+rest
    counts = defaultdict(int)
    for pairs in recurse(tuple(range(len(word)))):
        crossings = sum(a<c<b<d or c<a<d<b
                        for i, (a,b) in enumerate(pairs)
                        for c,d in pairs[i+1:])
        counts[crossings] += 1
    return dict(counts)


def moment(word, q, C):
    if len(word) % 2:
        return F(0)
    return C**(len(word)//2)*sum((count*q**k for k,count in
                                colored_pairings(tuple(word)).items()), F(0))


def integrate_poly(poly, q, C):
    return sum((coef*moment(('h',)*n,q,C) for n,coef in enumerate(poly)),F(0))


def block(abc, color='u'):
    a,b,c=abc
    return (color,)+('h',)*a+(color,)+('h',)*b+(color,)+('h',)*c+(color,)


def ladder_image(poly, basis, q, C, which):
    coeffs = basis_coefficients(poly,basis)
    out = [F(0)]
    for n, coeff in enumerate(coeffs):
        if which == 'P':
            out = add(out,basis[n],coeff*q**n)
        elif which == 'N':
            out = add(out,basis[n],coeff*n)
        elif which == 'D2' and n>=2:
            out = add(out,basis[n-2],coeff*sum(q**j for j in range(n))
                      *sum(q**j for j in range(n-1)))
    return out


def check(q,C):
    basis = hermites(q,C,12)
    configs = [abc for abc in product(range(7),repeat=3) if sum(abc)<=6]
    polynomials = {abc:conditional(block(abc),q,C,basis) for abc in configs}
    trace_checks = 0
    for abc, poly in polynomials.items():
        for n in range(5):
            direct = moment(('h',)*n+block(abc),q,C)
            expected = integrate_poly([F(0)]*n+poly,q,C)
            assert direct == expected,(q,C,abc,n)
            trace_checks += 1
    # The same h is retained in both blocks, whereas all u / v contractions
    # stay within their own four-endpoint block. Include h gaps between them.
    short = [abc for abc in configs if sum(abc)<=3]
    replica_checks = 0
    for abc in short:
        for xyz in short:
            for gap in range(3):
                direct = moment(tuple(reversed(block(abc)))+('h',)*gap+block(xyz,'v'),q,C)
                expected = integrate_poly([F(0)]*gap+mul(polynomials[abc],polynomials[xyz]),q,C)
                assert direct == expected,(q,C,abc,xyz,gap)
                replica_checks += 1
    # Independently act with W_n(alpha h+beta u), retaining the two outer u's.
    alpha,beta=F(3,5),F(4,5)
    ou_checks = 0
    for n in range(11):
        state = {('u',): F(1)}
        projected = [F(0)]
        for degree,coeff in enumerate(basis[n]):
            for tensor,value in fock_step(state,'u',q,C).items():
                if all(letter=='h' for letter in tensor):
                    projected=add(projected,basis[len(tensor)],coeff*value)
            nextstate=defaultdict(F)
            for color,scale in [('h',alpha),('u',beta)]:
                for tensor,value in fock_step(state,color,q,C).items():
                    nextstate[tensor]+=scale*value
            state={k:v for k,v in nextstate.items() if v}
        expected=add([F(0)],basis[n],C*q**n*alpha**n)
        if n>=2:
            expected=add(expected,basis[n-2],beta**2*C*C*
                         sum(q**j for j in range(n))*sum(q**j for j in range(n-1))*alpha**(n-2))
        assert not any(add(projected,expected,-1)),(q,C,n)
        ou_checks+=1
    # I_s = sum_n i^n s^(n+2)/(n+2)! sum_{a+b+c=n} Q_abc(h).
    # The proposed closed iota expression instead acts on x^(n+2), with i^(n+2).
    taylor=[]
    for n in range(7):
        total=[F(0)]
        for abc in configs:
            if sum(abc)==n:
                total=add(total,polynomials[abc])
        power=[F(0)]*(n+2)+[F(1)]
        diff=add(ladder_image(power,basis,q,C,'N'),power,-(n+2))
        proxy=add(ladder_image(diff,basis,q,C,'P'),[F(0)])
        proxy=[C*v/2 for v in proxy]
        proxy=add(proxy,ladder_image(power,basis,q,C,'D2'),-C*C)
        assert not any(add(total,proxy)),(q,C,n)
        taylor.append(dict(H_degree=n,coefficient_polynomial=list(map(str,total))))
    return dict(q=str(q),variance=str(C),conditional_words=len(configs),
                conditional_trace_checks=trace_checks,three_color_replica_checks=replica_checks,
                OU_Wick_degrees=ou_checks,Taylor_checks=len(taylor),Taylor_polynomials=taylor)


if __name__=='__main__':
    start=time.time()
    OUT.mkdir(parents=True,exist_ok=True)
    dest=OUT/'four-copy-q-checks.json'
    if dest.exists():raise FileExistsError(dest)
    rows=[]
    for q,C in product([F(1,4),F(1,2),F(3,4)],[F(1),F(3,2)]):
        rows.append(check(q,C))
        print(json.dumps(dict(q=str(q),C=str(C),status='passed')),flush=True)
    result=dict(rows=rows,elapsed_seconds=time.time()-start,
                created_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                environment=dict(python=platform.python_version(),arithmetic='fractions.Fraction'),
                scope='Exact q-Fock conditional polynomials versus independent colored pairings, three-color block replica factorization, same-copy OU Wick identity and four-u Taylor background. No finite-N fit, asymptotic proof or independent-probe replacement of the real charged force.',
                input_sha256={'code/memory_four_copy_background.py':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    with dest.open('x',encoding='utf-8') as stream:stream.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status='passed',parameter_sets=len(rows),
        checks={key:sum(row[key] for row in rows) for key in
        ['conditional_trace_checks','three_color_replica_checks','OU_Wick_degrees','Taylor_checks']},
        elapsed_seconds=result['elapsed_seconds'])))
