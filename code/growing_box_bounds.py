"""Finite-point checks of growing-box proof ingredients, not supremum certification."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,sys,platform
import mpmath

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/stage3/growing-box'

def ratio_log(c,Z,a,b):
    return c.loggamma(Z+a)-c.loggamma(Z+b)+(b-a)*c.log(Z)

def gauss_derivative_bound(c,A,B,C,z,j):
    assert c.re(B)>0 and c.re(C-B)>0 and 0<=z<1
    actual=c.rf(A,j)*c.rf(B,j)/c.rf(C,j)*c.hyp2f1(A+j,B+j,C+j,z)
    bound=(abs(c.rf(A,j)*c.gamma(C)/c.gamma(B)/c.gamma(C-B))
        *c.beta(c.re(B)+j,c.re(C-B))*(1-z)**(-max(c.re(A)+j,0)))
    return actual,bound

def main():
    OUT.mkdir(parents=True,exist_ok=True);dest=OUT/'diagnostics.json'
    if dest.exists():raise FileExistsError(dest)
    c=mpmath.mp.clone();c.dps=80
    checks=[];ratios=[];eulers=[];energies=[];tails=[]
    def check(name,error,tolerance):
        checks.append(dict(name=name,error=str(error),tolerance=str(tolerance),passed=bool(error<c.mpf(tolerance))))
    # The normalized ratio fixture bound is empirical; the theorem has an unspecified universal constant.
    for Q in [1,4,16,64]:
        for sign in [-1,1]:
            for factor in [16,64]:
                M=c.mpf(factor)*(1+Q)**2;Z=sign*2*c.j*M
                a=c.mpf('.5')+c.j*c.mpf('.7')*Q;b=c.mpf('.75')-c.j*c.mpf('.4')*Q
                low=ratio_log(c,Z,a,b)
                c.dps=110;high=ratio_log(c,Z,a,b);c.dps=80
                err=abs(low-high);scaled=abs(high)*M/(1+Q)**2
                check(f'gamma-precision-{Q}-{sign}-{factor}',err,'1e-65')
                check(f'gamma-fixture-scaled-bound-{Q}-{sign}-{factor}',scaled,'10')
                ratios.append(dict(Q=Q,M=str(M),sign=sign,log_ratio=[str(c.re(high)),str(c.im(high))],
                    relative_error=str(abs(c.expm1(high))),scaled_log_error=str(scaled),precision_change=str(err)))
    for Q in [1,4,16]:
        for tube in ['-.02','.02']:
            b=c.mpf('.5')+c.j*c.mpf('.6')*Q;d=c.mpf('.5')-c.j*c.mpf('.3')*Q
            delta=c.mpf(tube)+c.j*c.mpf('.2')*Q;g=b+d+delta
            triples=[('original',b,d,g),('plus_connection',b,d,1-delta),('minus_connection',g-b,g-d,1+delta)]
            for name,A,B,C in triples:
                for ztext in ['.2','.6']:
                    z=c.mpf(ztext)
                    for j in range(6):
                        actual,bound=gauss_derivative_bound(c,A,B,C,z,j)
                        ratio=abs(actual)/bound
                        check(f'euler-bound-{Q}-{tube}-{name}-{ztext}-{j}',ratio,'1.000000000000000000000000000001')
                        eulers.append(dict(Q=Q,tube_real_delta=tube,branch=name,z=ztext,derivative=j,
                            value=[str(c.re(actual)),str(c.im(actual))],bound=str(bound),ratio=str(ratio)))
    # Check the common inverse-mass shift cancellation using exact high-precision energies.
    for Mtext in ['100','10000','1000000']:
        M=c.mpf(Mtext);ell=c.log(2*M);L=3*c.log(ell);R=c.sqrt(24*c.log(ell));Q=L+R
        r=L*c.mpf('.7');x=R*c.mpf('.6');y=-R*c.mpf('.4');z=R*c.mpf('.2')
        energies0=[M+r+x,M+r+y,M+r+z,M+r]
        labels=[c.sqrt(v*v-c.mpf('.25')) for v in energies0]
        ch_m=x-y+z;ch_s=labels[0]-labels[1]+labels[2]-labels[3]
        common=[labels[i]-(energies0[i]-1/(8*M)) for i in range(4)]
        # Direct analytic bound: |s(m)-m+1/(8M)| <= Q/[4M^2]+1/[M^3] when |m-M|<=Q and M>2Q+1.
        assert M>2*Q+1
        bound=Q/(4*M*M)+1/(M**3)
        check(f'energy-offset-bound-{Mtext}',max(abs(v) for v in common)/bound,'1.00000000000000000001')
        check(f'energy-phase-difference-bound-{Mtext}',abs(ch_s-ch_m)/(4*bound),'1.00000000000000000001')
        energies.append(dict(M=Mtext,L=str(L),R=str(R),chi_mass=str(ch_m),chi_label=str(ch_s),
            difference=str(ch_s-ch_m),common_shift_remainders=[str(v) for v in common],each_remainder_bound=str(bound)))
    # Exact Gamma survival versus the exponential-moment Markov bound at the prescribed cutoffs.
    for q in [3,8,20]:
        for ell in [8,32,128]:
            beta0=c.mpf('.5');sigma=c.one
            L=(q+2)/beta0*c.log(ell);R=c.sqrt(8*(q+2)*c.log(ell))/sigma
            gaussian=c.erfc(sigma*R/2)**c.mpf('.5')
            gaussian_bound=c.exp(-sigma*sigma*R*R/8)
            check(f'Gaussian-tail-shape-{q}-{ell}',gaussian/gaussian_bound,'1.00000000000000000001')
            for p in [0,1,2]:
                exact=c.gammainc(p+1,L,c.inf)/c.factorial(p)
                markov=2**(p+1)*c.exp(-beta0*L)
                check(f'Gamma-moment-tail-{p}-{q}-{ell}',exact/markov,'1.00000000000000000001')
                tails.append(dict(q=q,ell=ell,p=p,L=str(L),R=str(R),Gamma_survival=str(exact),Markov_bound=str(markov),
                    normalized_gaussian_sqrt_tail=str(gaussian),normalized_gaussian_bound=str(gaussian_bound)))
    paths=['code/growing_box_bounds.py','sources/dlmf-5.11.html','sources/dlmf-15.8.html',
           'sources/kolchmeyer-liu-2411.08090v1.html','manuscript/endpoint.tex','manuscript/assembly.tex','manuscript/twobranch.tex']
    result=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),status='passed' if all(x['passed'] for x in checks) else 'failed',
        checks=checks,gamma_ratio_rows=ratios,Euler_derivative_rows=eulers,energy_rows=energies,cutoff_tail_rows=tails,
        parameters=dict(gamma_dps=[80,110],Euler_dps=80,nu='fixed by analytic argument; Euler fixtures are admissible parameter families',
            gamma_Q=[1,4,16,64],Euler_Q=[1,4,16],tube_real_delta=['-.02','.02'],Euler_derivatives=list(range(6)),
            cutoff_q=[3,8,20],cutoff_ell=[8,32,128],beta0='.5',sigma='1'),
        versions=dict(python=sys.version,mpmath=mpmath.__version__,platform=platform.platform()),
        limitations=['finite-point proof-ingredient checks, not uniform suprema or independent proof certification',
            'no exact finite-M 6j packet integral','no numerical certification of constants c,C,C0,kappa or minimum usable mass',
            'normalized Gamma fixture threshold 10 is only an empirical diagnostic, not the theorem constant',
            'the Gaussian tail check factors out the unknown positive-kernel envelope normalization'],
        hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths})
    dest.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(result['status'],len(checks),'checks',flush=True)
    print('max Gamma scaled log error',max(float(r['scaled_log_error']) for r in ratios),flush=True)
    print('max Euler ratio',max(float(r['ratio']) for r in eulers),flush=True)
    for x in checks:
        if not x['passed']:print('FAILED',x,flush=True)
    if result['status']!='passed':raise SystemExit(1)
if __name__=='__main__':main()
