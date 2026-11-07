"""Independent nonsymmetric four-arc sum at 50 decimal places."""
from pathlib import Path
import hashlib,json,platform
import mpmath as mp
from singlet_general_fourpoint import GeneralKernel
ROOT=Path(__file__).resolve().parent.parent

def reference(lam,times,L):
    q=mp.exp(-4*lam);r=q*q;T=mp.matrix(L)
    for n in range(L-1):
        T[n+1,n]=1/mp.sqrt(lam)
        T[n,n+1]=(1-q**(n+1))/(1-q)/mp.sqrt(lam)
    a,b,c,d=times;arcs=[d-c,c-b,b-a,a-d]
    U=[mp.expm(1j*s*T) for s in arcs];D=mp.diag([r**j for j in range(L)])
    f=[mp.mpf(1)];g=[mp.mpf(1)]
    for j in range(1,L):
        f.append(f[-1]*(1-q**j));g.append(g[-1]*(1-r*r*q**(j-1)))
    M=U[2]*D*U[1];terms=[]
    for i in range(L):
        for l in range(L-i):
            P=f[l+i]*g[i]/(f[l]*f[i])
            for m in range(L-i):
                terms.append(P*r**(i+l+m)*U[3][0,m+i]*M[m,l]*U[0][l+i,0])
    cross=q**4*mp.fsum(terms)
    unc=[]
    for s in [0,1]:
        j=[(i+s)%4 for i in range(4)]
        unc.append((U[j[2]]*D*U[j[1]]*U[j[3]]*D*U[j[0]])[0,0])
    return unc+[cross],mp.fsum(abs(v) for v in terms)/abs(mp.fsum(terms))

def main():
    target=ROOT/'results/susskind-first-order/representation-checks.json'
    if target.exists():raise FileExistsError(target)
    mp.mp.dps=50;rows=[]
    for ls in ['0.25','0.1']:
        for ts in [['.3','.1','-.1','0'],['.2','-.05','.35','.1']]:
            ref,cond=reference(mp.mpf(ls),list(map(mp.mpf,ts)),14)
            val=GeneralKernel(float(ls),14).four(list(map(float,ts)))
            errs=[abs(complex(z)-complex(re,im)) for z,re,im in zip(ref,val['channels'],val['channel_imag'])]
            assert max(errs)<1e-12
            rows.append(dict(lambda_value=ls,times=ts,L=14,reference=[str(z) for z in ref],absolute_errors=errs,condition=str(cond)))
    names=['code/check_general_fourpoint.py','code/singlet_general_fourpoint.py','code/singlet_otoc.py']
    data=dict(status='passed',scope='Same finite L only; independent nonsymmetric BINT index sum versus symmetric Jacobi bilinear channels',dps=50,records=rows,
        environment=dict(python=platform.python_version(),mpmath=mp.__version__),input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in names})
    target.write_text(json.dumps(data,indent=2)+'\n');print('Passed 12 channel comparisons at 50 digits.')
if __name__=='__main__':main()
