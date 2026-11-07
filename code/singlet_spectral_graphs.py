"""Fixed-order centered singlet graphs; all allowed Gaussian pairings retained."""
from itertools import combinations
from math import comb
from singlet_coefficient import pairings, dual_product


def minus(a, b):
    return tuple(x-y for x, y in zip(a, b))


def graph(length, pairing, marks):
    pi = {a:b for a,b in pairing} | {b:a for a,b in pairing}
    def entry(a,b):
        return (1 if a <= pi[b] else -1,
                -2*int((pi[b], a) == marks[0]),
                -2*int((pi[b], a) == marks[1]), 0)
    singles = [entry(i,i) for i in range(length)]
    assert all(s[1:] == (0,0,0) for s in singles)
    kappa = [-comb(length,2)-length//2, 0, 0, 0]
    for i,j in combinations(range(length), 2):
        minor = minus(dual_product(entry(i,i),entry(j,j)),
                      dual_product(entry(i,j),entry(j,i)))
        den = singles[i][0]*singles[j][0]
        kappa = [a+b*den for a,b in zip(kappa,minor)]
    cr = sum(a<c<b<d or c<a<d<b for i,(a,b) in enumerate(pairing) for c,d in pairing[i+1:])
    assert kappa[0] == -4*cr
    return dict(pairing=pairing, marks=marks, crossings=cr,
                alpha=kappa[1], beta=kappa[2], delta=kappa[3])


def centered_graphs(order):
    """tr O (ad_H)^order O; order must be even, including zero."""
    assert order % 2 == 0
    length = order+4
    rows=[]
    for k in range(order+1):
        marks=((0,1),(order-k+2,order-k+3))
        for pairing in pairings(tuple(range(length))):
            # Wick centering at the two A insertions removes exactly these.
            if any(m in pairing for m in marks):
                continue
            row=graph(length,pairing,marks)
            row.update(k=k, sign=(-1)**k*comb(order,k))
            rows.append(row)
    return rows


def aggregate(rows):
    out={}
    for r in rows:
        c=r['crossings']
        values=out.setdefault(c,[0,0])
        values[0]+=r['sign']*r['delta']
        values[1]+=r['sign']*r['alpha']*r['beta']
    return out


if __name__ == '__main__':
    for order in [0,2,4,6]:
        rows=centered_graphs(order)
        print(order,len(rows),aggregate(rows),flush=True)
