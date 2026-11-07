"""Exact charged covariance lift, marked products, and L2 rank-filter bridge."""
from pathlib import Path
from itertools import product
from math import comb
import hashlib,json,platform
import numpy as np
import mpmath as mp
from singlet_wick_moments import gue_basis
from cubic_covariance_bridge import ratio
from cubic_linear_limit import distribution,scaled_coefficient
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-filtered-spectrum'

def finite(N,p,m):
    cs,v,r=gue_basis(N,p);_,u,s=gue_basis(N,p-1);dim=1<<N
    charges=np.array([x.bit_count() for x in range(dim)]);Q=np.diag(charges)
    E=[[ci.conj().T@cj for cj in cs] for ci in cs]
    def cas2(a):
        out=np.zeros_like(a)
        for i in range(N):
            for j in range(N):
                b=E[j][i]@a-a@E[j][i]
                out+=E[i][j]@b-b@E[i][j]
        return out
    ix=np.flatnonzero(charges==m);vm=v[:,ix[:,None],ix];um=u[:,ix[:,None],ix]
    numer=(N-m+p)*(2*p*(N-p+1)*v-np.array([cas2(a) for a in v]))
    kn=(m-p+1)*p*p
    fm=numer[:,ix[:,None],ix]
    Cp=comb(m,p)*comb(N-m+p,p);Cr=comb(m,p-1)*comb(N-m+p-1,p-1)
    channels=[]
    for n in [m,m-1,m-2]:
        if n<p-1:continue
        iy=np.flatnonzero(charges==n);vn=v[:,iy[:,None],iy];un=u[:,iy[:,None],iy]
        left=np.einsum('a,aij,akl->ijkl',r,fm,vn)
        right=np.einsum('a,aij,akl->ijkl',s,um,un)
        # Left denominator 4 kn Cp; right denominator 2 Cr and lift (n-p+1)/(m-p+1).
        l=left*Cr*(m-p+1);rr=right*(2*kn*Cp*(n-p+1))
        assert max(float(np.max(abs(l))),float(np.max(abs(rr))))<2**51
        assert np.array_equal(l,rr),(N,p,m,n)
        channels.append(dict(output_charge=n,entries=int(l.size)))
    D=lambda a:Q@a-sum(c.conj().T@a@c for c in cs)
    h=sum(((i%3)-1)*a for i,a in enumerate(v));g=sum(((i%5)-2)*a for i,a in enumerate(u))
    # Integer numerator of the one-step p-body lift g -> (Q-p+1)g/(m-p+1).
    lift_num=(Q-(p-1)*np.eye(dim))@g
    checks=0
    for degree in [2,3]:
        for at in range(degree):
            arr=[h]*degree;arr[at]=g
            prod0=np.eye(dim,dtype=complex);prod1=prod0.copy()
            for j,a in enumerate(arr):prod0=prod0@a;prod1=prod1@(lift_num if j==at else a)
            lhs=(degree*p*prod1-D(prod1))[np.ix_(ix,ix)]
            natural=((degree*p-1)*prod0-D(prod0))[np.ix_(ix,ix)]
            rhs=(m-p)*natural+(degree-1)*p*prod0[np.ix_(ix,ix)]
            assert np.max(abs(lhs))<2**45 and np.array_equal(lhs,rhs)
            checks+=1
    return dict(N=N,p=p,m=m,covariance_blocks=channels,marked_product_checks=checks)

def error(lam,p,dps):
    mp.mp.dps=dps;lam=mp.mpf(lam);N=int(p*p/lam);m=N//2;w=distribution(N,p,m)
    co=[scaled_coefficient(N,p,m,p-s) for s in range(p+1)]
    avg=mp.fsum(x*y for x,y in zip(w,co));ff=[x-avg for x in co]
    var=mp.fsum(x*y*y for x,y in zip(w,ff));ff=[x/mp.sqrt(var) for x in ff]
    rr=[ratio(N,p,m,s,1) for s in range(p+1)]
    yy=[mp.sqrt(lam)*(mp.mpf(x.numerator)/x.denominator-1) for x in rr]
    err=mp.fsum(x*(y-z)**2 for x,y,z in zip(w,ff,yy))
    # Exact first energy-sum F-Y overlap, using the preceding finite Jordan coefficients.
    kap=[];sc=[];cp=mp.mpf(m-p)/m*mp.mpf(N-m+p+1)/(N-m+1)
    a0=2*p-m+m*cp
    for a in range(p+1):
        k=mp.fprod((1-mp.mpf(a)/(m-j))*(1-mp.mpf(a)/(N-m+p-j)) for j in range(p))
        jj=mp.mpf(p)*(p*(N-m+1)-a*(N-a+1))/((m-a)*(N-m-a+1))
        kap.append(k);sc.append(a0+mp.mpf(p*p)/m+k*jj)
    wa=list(reversed(w));mean=lambda v:mp.fsum(x*y for x,y in zip(wa,v))
    cn=(mean(sc)-a0)/(1+mean(kap));gv=[s-cn*k for s,k in zip(sc,kap)];bg=mean(gv)
    # Original unnormalized dimensionless F coefficients reverse ff and restore sqrt(var).
    fg=mp.fsum(wa[a]*ff[p-a]*mp.sqrt(var)*(gv[a]-bg) for a in range(p+1))/var
    c=2*lam*(1-mp.exp(-4*lam))/(1+mp.exp(-4*lam))
    return dict(N=N,p=p,lambda_value=str(lam),precision=dps,L2_bridge_error_squared=mp.nstr(err,60),
        p_F_JY_over_Fnorm=mp.nstr(p*fg,60),limit_p_F_JY_over_Fnorm=mp.nstr((1-2*lam-c)/2,60))

def main():
    OUT.mkdir(parents=True,exist_ok=True);target=OUT/'charged-bridge.json'
    if target.exists():raise FileExistsError(target)
    finite_rows=[finite(*x) for x in [(4,2,2),(5,2,3),(6,2,3)]]
    rows=[];changes=[]
    for lam in ['0.1','0.25','0.5']:
        for p in [16,32,64,128]:
            a=error(lam,p,80);b=error(lam,p,120)
            change=max(abs(mp.mpf(a[k])-mp.mpf(b[k])) for k in ['L2_bridge_error_squared','p_F_JY_over_Fnorm'])
            assert change<mp.mpf('1e-55');changes.append(change);rows.append(b)
            print(lam,p,b['L2_bridge_error_squared'][:12],b['p_F_JY_over_Fnorm'][:14],flush=True)
    names=['code/singlet_filtered_bridge.py','code/cubic_covariance_bridge.py','code/cubic_linear_limit.py',
        'code/cubic_linear_dynamics.py','code/cubic_linear_channels.py','code/singlet_wick_moments.py','code/complex_syk_dictionary.py']
    target.write_text(json.dumps(dict(status='passed',finite=finite_rows,rank_records=rows,
        max_precision_change=mp.nstr(max(changes),25),precisions=[80,120],
        scope='Exact bounded Gaussian-integer charged tensors and marked identities; full-rank L2 sums; empirical precision checks',
        environment=dict(python=platform.python_version(),numpy=np.__version__,mpmath=mp.__version__),
        input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names}),indent=2)+'\n')
if __name__=='__main__':main()
