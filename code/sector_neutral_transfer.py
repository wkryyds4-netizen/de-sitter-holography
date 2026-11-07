"""Pure finite Hermite/filter and selected local-descendant interfaces."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import comb, factorial, sqrt
import numpy as np
from complex_syk_dictionary import fermions, creators
import memory_wick_descendant_chords as dc
import slater_higher_wick_sources as src


def mu(n, p, m, rank):
    choose = lambda a, b: comb(a, b) if 0 <= b <= a else 0
    return F(choose(m-rank, p)*choose(n-m+p-rank, p),
             comb(m, p)*comb(n-m+p, p))


def finite_system(n, p, m):
    cs = fermions(n)
    creation = creators(cs, p)
    ids = [i for i in range(1 << n) if i.bit_count() == m]
    variance = comb(m, p)*comb(n-m+p, p)
    basis = []
    for a, ca in enumerate(creation):
        basis.append((ca@ca.conj().T)[np.ix_(ids, ids)]/sqrt(variance))
        for cb in creation[:a]:
            x = (ca@cb.conj().T)[np.ix_(ids, ids)]
            basis.extend([(x+x.conj().T)/sqrt(2*variance),
                          1j*(x-x.conj().T)/sqrt(2*variance)])
    d = len(ids)
    eye = np.eye(d)
    casimir = np.zeros((d*d, d*d), complex)
    E = [[(a.conj().T@b)[np.ix_(ids, ids)] for b in cs] for a in cs]
    for i in range(n):
        for j in range(n):
            left = np.kron(E[i][j], eye)-np.kron(eye, E[i][j].T)
            right = np.kron(E[j][i], eye)-np.kron(eye, E[j][i].T)
            casimir += left@right
    ev, vectors = np.linalg.eigh(casimir)
    ranks = np.array([min(range(min(m,n-m)+1),
                         key=lambda r: abs(x-2*r*(n+1-r))) for x in ev])
    assert max(abs(ev-2*ranks*(n+1-ranks))) < 1e-10
    phi = sum(np.kron(b, b.T) for b in basis)
    spectral = (vectors*np.array([float(mu(n,p,m,int(r))) for r in ranks]))@vectors.conj().T
    return basis, vectors, ranks, phi, float(np.linalg.norm(phi-spectral))


def hermite_add(a, b, factor=1):
    out = {k: x.copy() for k, x in a.items()}
    for k, x in b.items():
        out[k] = out.get(k, 0)+factor*x
    return out


def hermite_step(poly, basis, creation_only=False):
    out = {}
    for key, x in poly.items():
        for i, b in enumerate(basis):
            y = b@x
            up = list(key)
            up[i] += 1
            up = tuple(up)
            out[up] = out.get(up, 0)+y
            if key[i] and not creation_only:
                down = list(key)
                down[i] -= 1
                down = tuple(down)
                out[down] = out.get(down, 0)+key[i]*y
    return out


def filter_rows(n, p, m, maximum):
    basis, vectors, ranks, phi, spectral_error = finite_system(n,p,m)
    d = len(basis[0])
    zero = (0,)*len(basis)
    Z = {zero: np.eye(d, dtype=complex)}
    Q, previous = Z, {}
    lipschitz = F(p, m-p+1)+F(p,n-m+1)
    rows = []
    for j in range(maximum+1):
        E = hermite_add(Z,Q,-1)
        e2 = de2 = dde2 = filter2 = direct_error2 = 0.
        chaos0 = 0.
        for key, x in E.items():
            degree = sum(key)
            weight = float(np.prod([factorial(k) for k in key]))/d
            raw = x.reshape(-1)
            spectral = vectors.conj().T@raw
            deficit = degree*p-ranks
            multiplier = np.array([float(mu(n,p,m,int(r))-mu(n,p,m,degree*p))
                                   for r in ranks])
            eta_filter = p*multiplier*spectral
            direct = p*(phi@raw-float(mu(n,p,m,degree*p))*raw)
            e2 += weight*float(np.vdot(raw,raw).real)
            de2 += weight*float(np.vdot(deficit*spectral,deficit*spectral).real)
            dde2 += weight*float(np.vdot(deficit**2*spectral,deficit**2*spectral).real)
            filter2 += weight*float(np.vdot(eta_filter,eta_filter).real)
            direct_error2 += weight*float(np.linalg.norm(direct-vectors@eta_filter)**2)
            if degree == 0:
                chaos0 += weight*float(np.vdot(raw,raw).real)
        interpolation_rhs = sqrt(e2*dde2)
        filter_rhs = float(p*lipschitz)**2*de2
        assert de2 <= interpolation_rhs+1e-10*max(1,interpolation_rhs)
        assert filter2 <= filter_rhs+1e-10*max(1,filter_rhs)
        assert direct_error2 < 1e-18
        rows.append(dict(degree=j, hermite_matrices=len(E), E_norm_squared=e2,
                         deficit_E_norm_squared=de2, second_deficit_E_norm_squared=dde2,
                         interpolation_rhs=interpolation_rhs,
                         eta_filter_norm_squared=filter2, filter_bound_squared=filter_rhs,
                         direct_spectral_error=sqrt(direct_error2),
                         scalar_chaos_E_norm_squared=chaos0))
        if j < maximum:
            an=sum(float(mu(n,p,m,v*p)) for v in range(j))
            Qnext=hermite_add(hermite_step(Q,basis),previous,-an)
            previous,Q=Q,Qnext
            Z=hermite_step(Z,basis,True)
    return dict(N=n,p=p,m=m,C=1,Gaussian_coordinates=len(basis),maximum=maximum,
                phi_casimir_error=spectral_error,lipschitz=str(lipschitz),rows=rows)


def peval(poly, gamma, lam):
    return sum((v*lam**a*gamma**(4*b) for (a,b),v in poly.items()), F(0))


@lru_cache(None)
def term_word(terms, gamma, lam, mode):
    choices = []
    for term in terms:
        degrees = [term.selected_chaos] if term.selected_chaos >= 0 else range(term.length%2,term.length+1,2)
        choices.append(degrees if term.ou_power or term.selected_chaos >= 0 else [None])
    total=F(0)
    for degrees in product(*choices):
        offset=0
        marks=[];wicks=[];chaoses=[]
        factor=F(1)
        for term, degree in zip(terms,degrees):
            marks += [(a+offset,b+offset,'L') for a,b in term.mark_blocks]
            wicks += [(a+offset,b+offset) for a,b in term.wick_blocks]
            if degree is not None:
                if term.length:
                    chaoses.append((offset,offset+term.length,degree))
                elif degree != 0:
                    factor=F(0)
                factor *= gamma**(4*term.ou_power*degree)
            offset += term.length
        if factor:
            total += factor*src.moment(offset,marks,gamma,lam,wicks,chaoses,mode=mode)
    return total


def field_word(fields,gamma,lam,mode='trace',dagger_first=False):
    gamma,lam=F(gamma),F(lam)
    total=F(0)
    for selected in product(*(list(field.terms.items()) for field in fields)):
        terms=[item[0] for item in selected]
        if dagger_first:
            terms[0]=terms[0].dagger()
        coefficient=F(1)
        for _,poly in selected:
            coefficient *= peval(poly,gamma,lam)
        total += coefficient*term_word(tuple(terms),gamma,lam,mode)
    return total


def poly_add(a,b,factor=F(1)):
    out=dict(a)
    for k,v in b.items():
        out[k]=out.get(k,F(0))+factor*v
    return {k:v for k,v in out.items() if v}


def xmul(poly):
    return {k+1:v for k,v in poly.items()}


def regression_fields(maximum,gamma,lam):
    """C=1 regressions from the analytic forced q-Hermite recurrence."""
    q=F(gamma)**4;lam=F(lam)
    qint=lambda j: sum((q**k for k in range(j)),F(0))
    H=[{0:F(1)},{1:F(1)}]
    for j in range(1,maximum):
        H.append(poly_add(xmul(H[j]),H[j-1],-qint(j)))
    kappas=[F(0)]
    derivative=F(0)
    for j in range(1,maximum+1):
        derivative += sum((k*q**k for k in range(j)),F(0))/qint(j)
        kappas.append(lam*(j*j-4*derivative))
    g=[{} for _ in range(maximum+1)]
    for j in range(2,maximum):
        forcing=-4*lam*sum((q**v*kappas[v] for v in range(j)),F(0))
        g[j+1]=poly_add(poly_add(xmul(g[j]),g[j-1],-qint(j)),H[j-1],forcing)
    vfields=[]
    for j in range(maximum+1):
        poly={k:4*lam*q**j*kappas[j]*v for k,v in H[j].items()}
        rest=dict(g[j]);expansion={}
        for degree in reversed(range(j+1)):
            c=rest.get(degree,F(0))
            if c:
                expansion[degree]=c
                rest=poly_add(rest,H[degree],-c)
        assert not rest
        for degree,c in expansion.items():
            poly=poly_add(poly,H[degree],(q**j-q**degree)*c)
        vfields.append(dc.Field(j,{dc.Term(k):dc.mono(coefficient=c) for k,c in poly.items()},f'v{j}').validate())
    return vfields,kappas,g
