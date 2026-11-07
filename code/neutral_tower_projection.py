"""Pure one-U projection calculus, C=1; no file writes or scientific mains."""
from fractions import Fraction as F
from functools import lru_cache
from neutral_selected_response import hermite,pair_moment,power_coeff

def add(a,b,scale=1):
    out=dict(a)
    for k,v in b.items():out[k]=out.get(k,F(0))+scale*v
    return {k:v for k,v in out.items() if v}

def qi(n,q):return sum((q**i for i in range(n)),F(0))

def kappa(n,q,lam,charged=False):
    d=sum((sum((j*q**j for j in range(k)),F(0))/qi(k,q) for k in range(1,n+1)),F(0))
    return lam*(n*(n-int(charged))-4*d)

def left(module,q):
    """h times Wick(h^a U h^b), in the one-U Wick basis."""
    out={}
    for (a,b),v in module.items():
        out=add(out,{(a+1,b):v})
        if a:out=add(out,{(a-1,b):v*qi(a,q)})
        if b:out=add(out,{(a,b-1):v*q**(a+1)*qi(b,q)})
    return out

@lru_cache(None)
def symbol(a,b,q):
    """Wick(h^a U h^b)=P_ab(h_L,h_R)U; q-Hermite parameter r=q."""
    if a==0:return {(0,j):v for j,v in hermite(b,q).items()}
    out={(i+1,j):v for (i,j),v in symbol(a-1,b,q).items()}
    if a>1:out=add(out,symbol(a-2,b,q),-qi(a-1,q))
    if b:out=add(out,symbol(a-1,b-1,q),-q**a*qi(b,q))
    return out

def to_symbol(module,q):
    out={}
    for (a,b),v in module.items():out=add(out,symbol(a,b,q),v)
    return out

def tower(maximum,q,lam):
    """P eta_j=g_j+sqrt(lambda) E_j; P N_j=sqrt(lambda) Xi_j."""
    E=[{} for _ in range(maximum+1)];g=[{} for _ in E]
    for j in range(2,maximum):
        source={};reg=F(0)
        for v in range(1,j):
            source=add(source,{(i,j-2-i):F(1) for i in range(v)},4*lam*q**v)
            reg+=q**v*kappa(v,q,lam)
        E[j+1]=add(add(left(E[j],q),E[j-1],-qi(j,q)),source,-1)
        gg={}
        for d,c in g[j].items():
            gg=add(gg,{d+1:c})
            if d:gg=add(gg,{d-1:c*qi(d,q)})
        g[j+1]=add(add(gg,g[j-1],-qi(j,q)),{j-1:4*lam*reg},-1)
    Xi=[{}];vreg=[{}]
    for j in range(1,maximum+1):
        Xi.append(add({(i,j-1-i):4*lam*q**j for i in range(j)},
                      {(a,b):(q**j-q**(a+b+1))*v for (a,b),v in E[j].items()}))
        vreg.append(add({j:4*lam*q**j*kappa(j,q,lam)},
                        {d:(q**j-q**d)*v for d,v in g[j].items()}))
    return E,Xi,g,vreg

def energy_polynomial(reg,q):
    out={}
    for d,v in reg.items():out=add(out,hermite(d,q),v)
    return out

def pair_integral(poly,q,a=0,b=0):
    return sum((v*pair_moment(i+a,j+b,q) for (i,j),v in poly.items()),F(0))

def selected_symbol(poly,q,lam,gamma):
    """sigma(P(h_L,h_R) W), sqrt(q)=gamma^2 in rational tests."""
    return -gamma**2*lam*(pair_integral(poly,q,1,0)+pair_integral(poly,q,0,1))
