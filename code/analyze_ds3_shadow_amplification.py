"""Independent elementary modulus check for the two complete Omega coefficients."""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from run_ds3_collision import clean


def main():
    root=Path(__file__).resolve().parents[1];base=root/'results/stage4/ds3-full-tail'
    target=base/'amplification.json'
    if target.exists():raise RuntimeError('Preserve prior result')
    study=json.loads((base/'study.json').read_text());mp.mp.dps=150;rows=[]
    for row in study['rows']:
        if row['dps']!=150:continue
        M=mp.mpf(row['M']);u=mp.mpf('.3');v=mp.mpf('-.2');p=q=mp.mpf('.6')
        a=M+u+v;b=M+u;c=M;d=M+v
        xs=[(b+a-q)/2,(b-a+q)/2,(-b+p-c)/2,(-b+c-p)/2]
        common=mp.sinh(mp.pi*b)/(4*mp.pi**2*mp.sqrt(mp.fprod(mp.cosh(mp.pi*x) for x in xs)))
        exact=[common*mp.cosh(mp.pi*(sigma*d-a-p)/2)*mp.cosh(mp.pi*(sigma*d-q-c)/2) for sigma in [1,-1]]
        normal=mp.mpc(*row['details']['normalization'])
        saved=[abs(normal*mp.mpc(*sh['coefficient'])) for sh in row['details']['shadows']]
        relative=[abs(x-y)/x for x,y in zip(exact,saved)]
        asymptotic_ratio=mp.exp(2*mp.pi*M+mp.pi*(u+3*v+p+q)/2)/(4*mp.cosh(mp.pi*(u+p)/2)*mp.cosh(mp.pi*(v-q)/2))
        correction=(1+mp.exp(-mp.pi*(d+a+p)))*(1+mp.exp(-mp.pi*(d+q+c)))
        ratio=exact[1]/exact[0]
        plus_limit=mp.exp(mp.pi*(u-v+p+q)/4)*mp.sqrt(mp.cosh(mp.pi*(u+p)/2)*mp.cosh(mp.pi*(v-q)/2))/(4*mp.pi**2)
        x,y=row['details']['shadows'][1]['omega']
        # Error amplification when only Omega(h,1-hp) is perturbed; y is held at its stored sum.
        gain=exact[1]*abs(mp.mpc(*y['value']))
        rows.append(clean(dict(M=int(M),coefficient_moduli=exact,gamma_moduli=saved,relative_errors=relative,
                    shadow_ratio=ratio,asymptotic_ratio=asymptotic_ratio,exponential_correction=correction,
                    ratio_identity_error=abs(ratio/asymptotic_ratio-correction),plus_limit=plus_limit,
                    log10_shadow_ratio=mp.log10(ratio),gain_for_first_omega=gain,
                    first_omega_tail=mp.mpf(x['tail']),first_omega_tail_amplified=gain*mp.mpf(x['tail']),
                    first_omega_accuracy_for_1e_minus_6=mp.mpf('1e-6')/gain)))
    paths=['code/analyze_ds3_shadow_amplification.py','results/stage4/ds3-full-tail/study.json']
    result=dict(status='passed' if all(max(map(mp.mpf,r['relative_errors']))<mp.mpf('1e-140') and mp.mpf(r['ratio_identity_error'])<mp.mpf('1e-140') for r in rows) else 'failed',
                rows=rows,dps=150,hashes={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in paths})
    target.write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(result['status'])
    for r in rows:print(r['M'],'log10 ratio',mp.nstr(mp.mpf(r['log10_shadow_ratio']),9),'required Omega accuracy',mp.nstr(mp.mpf(r['first_omega_accuracy_for_1e_minus_6']),9))


if __name__=='__main__':main()
