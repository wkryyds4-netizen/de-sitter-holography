"""Positive Gauss rules for the normalized sin^4 spectral preparation."""
import mpmath
import numpy as np
from scipy.special import roots_legendre
from quantized_observer import mass


def packet_density(d,lo=.7,hi=.9):
    x=np.asarray(d,dtype=float);w=hi-lo
    return np.where((x>=lo)&(x<=hi),8/(3*w)*np.sin(np.pi*(x-lo)/w)**4,0.)


def gauss_packet(n,lo="0.7",hi="0.9",dps=70):
    """Stieltjes construction in x in [-1,1]; weight 4/3*cos(pi*x/2)^4.

    Returns full decimal strings. The first 2*n moments are matched to the
    requested arithmetic precision, not a certified bound for a general
    conditionally integrated function.
    """
    if not 1<=n<=12:raise ValueError("Supported weighted orders: 1..12")
    c=mpmath.mp.clone();c.dps=dps
    low,high=c.mpf(lo),c.mpf(hi)
    if low<=0 or high<=low:raise ValueError("Positive ordered spectral endpoints required")
    moments=[c.zero if k%2 else c.quad(lambda x,k=k:x**k*4*c.cos(c.pi*x/2)**4/3,[-1,0,1]) for k in range(2*n+1)]
    def inner(a,b):return c.fsum(ai*bj*moments[i+j] for i,ai in enumerate(a) for j,bj in enumerate(b))
    previous=[];poly=[c.one];previous_norm=c.one;off=[]
    for degree in range(n):
        norm=inner(poly,poly)
        beta=c.zero if degree==0 else norm/previous_norm
        if degree:off.append(c.sqrt(beta))
        # Symmetry makes every diagonal Jacobi coefficient zero exactly.
        following=[c.zero]+poly
        for i,x in enumerate(previous):following[i]-=beta*x
        previous,poly,previous_norm=poly,following,norm
    jac=c.matrix(n)
    for i,x in enumerate(off):jac[i,i+1]=jac[i+1,i]=x
    eigen,vectors=c.eigsy(jac)
    center=(low+high)/2;half=(high-low)/2
    nodes=[];weights=[]
    for i in range(n):
        x=c.zero if n%2 and i==n//2 else eigen[i]
        nodes.append(c.nstr(center+half*x,40))
        weights.append(c.nstr(moments[0]*vectors[0,i]**2,40))
    return dict(order=n,lo=str(lo),hi=str(hi),nodes=nodes,weights=weights,
        construction_dps=dps,polynomial_exact_degree=2*n-1,
        density="8/(3w)*sin(pi*(d-lo)/w)^4 on [lo,hi]",normalization="integral P(d) dd = 1")


def inner_rule(n,lo=.2,hi=3.):
    x,w=roots_legendre(n)
    return lo+(hi-lo)*(x+1)/2,w*(hi-lo)/2


def cartesian_tasks(n,rule,config,nu="0.6",lower=.2,upper=3.):
    from crossing_point_cache import point_identity,point_key
    nodes,weights=inner_rule(n,lower,upper)
    tasks=[]
    for di,d in enumerate(rule["nodes"]):
        for i,a in enumerate(nodes):
            for bi,b in enumerate(nodes):
                for ci in range(i,n):
                    z=nodes[ci]
                    labels=[format(float(v),".17g") for v in (a,b,z)]+[d]
                    identity=point_identity(config,nu,labels)
                    weight=(1 if i==ci else 2)*weights[i]*weights[bi]*weights[ci]
                    tasks.append(dict(key=point_key(identity),identity=identity,outer_index=di,inner_weight=float(weight)))
    return tasks


def rephase(tasks,records,rule,times=(0.,.5,1.,2.),sigma=1.):
    import math
    pieces=[[[] for t in times] for _ in rule["nodes"]]
    absolute=[0. for _ in rule["nodes"]]
    for task in tasks:
        record=records[task["key"]]
        raw=complex(float(record["conditional"]["real"]),float(record["conditional"]["imag"]))
        a,b,z,d=map(float,task["identity"]["labels"])
        ma,mb,mc,md=map(float,map(mass,(a,b,z,d)))
        gaussian=np.exp(-.5*sigma*sigma*((md-mc)**2+(mc-mb)**2+(mb-ma)**2+(ma-md)**2))
        base=raw*gaussian*task["inner_weight"]
        di=task["outer_index"];absolute[di]+=abs(base)
        for ti,t in enumerate(times):pieces[di][ti].append(base*np.exp(-1j*t*(md-mc+mb-ma)))
    conditional=[]
    for at_d in pieces:
        conditional.append([complex(math.fsum(z.real for z in xs),math.fsum(z.imag for z in xs)) for xs in at_d])
    weights=list(map(float,rule["weights"]))
    total=[complex(math.fsum(w*vals[i].real for w,vals in zip(weights,conditional)),
                   math.fsum(w*vals[i].imag for w,vals in zip(weights,conditional))) for i in range(len(times))]
    return total,conditional,absolute