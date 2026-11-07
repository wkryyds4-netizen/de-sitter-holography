"""High-precision nonsymmetric BINT transfer sum versus positive Jacobi form."""
from pathlib import Path
import hashlib,json,platform
import mpmath as mp
from singlet_otoc import Kernel

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-singlet-otoc'


def nonsymmetric(lam,t,L):
    q=mp.exp(-4*lam);r=q*q;T=mp.matrix(L)
    for n in range(L-1):
        T[n+1,n]=1/mp.sqrt(lam)
        T[n,n+1]=(1-q**(n+1))/(1-q)/mp.sqrt(lam)
    minus=mp.expm(-1j*t*T);plus=minus.conjugate()
    D=mp.diag([r**j for j in range(L)])
    M=minus*D*plus
    f=[mp.mpf(1)];g=[mp.mpf(1)]
    for j in range(1,L):f.append(f[-1]*(1-q**j));g.append(g[-1]*(1-r*r*q**(j-1)))
    terms=[]
    for i in range(L):
        for l in range(L-i):
            P=f[l+i]*g[i]/(f[l]*f[i])
            for m in range(L-i):
                terms.append(P*r**(i+l+m)*plus[0,m+i]*M[m,l]*minus[l+i,0])
    crossed=mp.fsum(terms)
    uncrossed=(minus*D*plus*plus*D*minus)[0,0]
    return crossed,uncrossed,mp.fsum(abs(v) for v in terms)/max(abs(crossed),mp.mpf('1e-100'))


def main():
    target=OUT/'representation-checks.json'
    if target.exists():raise FileExistsError(target)
    mp.mp.dps=50;rows=[]
    for lam in ['0.25','0.1']:
        kernel=Kernel(float(lam),18)
        for t in ['0.2','0.7']:
            c,u,condition=nonsymmetric(mp.mpf(lam),mp.mpf(t),18)
            v=kernel.at(float(t))
            dc=abs(float(c.real)-v['crossed_reduced']);du=abs(float(u.real)-v['uncrossed'])
            assert dc<1e-12 and du<1e-12 and abs(c.imag)<mp.mpf('1e-40') and abs(u.imag)<mp.mpf('1e-40')
            rows.append(dict(lambda_value=lam,time=t,L=18,crossed_mp=str(c),uncrossed_mp=str(u),
                crossed_absolute_difference=dc,uncrossed_absolute_difference=du,nonsymmetric_sum_condition=str(condition)))
    normalization=[]
    for z in ['0.025','0.05','0.1']:
        zz=mp.mpf(z);u=zz**-4*mp.hyperu(4,1,1/zz)
        integ=mp.quad(lambda y:y**3*mp.exp(-y)/(1+zz*y)**4,[0,1,mp.inf])/mp.gamma(4)
        assert abs(u-integ)<mp.mpf('1e-40')
        normalization.append(dict(z=z,tricomi=str(u),gamma_integral=str(integ),difference=str(abs(u-integ))))
    inputs=['code/check_otoc_representation.py','code/singlet_otoc.py']
    data=dict(status='passed',scope='Same finite L=18; independently verifies diagonal similarity, all three channel factors, and Tricomi normalization. Does not certify infinite-cutoff roundoff.',
        dps=50,records=rows,tricomi_normalization=normalization,environment=dict(python=platform.python_version(),mpmath=mp.__version__),
        input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs})
    target.write_text(json.dumps(data,indent=2)+'\n');print('4 high-precision channel pairs and 3 gamma integrals passed',flush=True)

if __name__=='__main__':main()
