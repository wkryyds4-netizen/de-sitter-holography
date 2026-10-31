"""Thermal limit and actual Gaussian-regulated vacuum-return diagnostics."""
import hashlib,json,platform
from pathlib import Path
import mpmath as mp
from ds3_gaussian_field import limiting,return_block,density,multiplier_bound
from ds3_fock_vacuum import block_bound
from run_ds3_collision import clean

def main():
    root=Path(__file__).resolve().parents[1];base=root/'results/stage4/ds3-fock-vacuum';dest=base/'gaussian.json';raw=base/'gaussian-runs.jsonl'
    if dest.exists() or raw.exists():raise RuntimeError('Preserve Gaussian evidence')
    mp.mp.dps=60;p=mp.mpf('.6');sigma=mp.mpf('.8');limits=[];rows=[];checks=[]
    for t in map(mp.mpf,['0','.5','1','2']):
        L=limiting(t,p,sigma);der=limiting(t,p,sigma,derivative=1)
        reversed_value=limiting(-t,p,sigma);continued=limiting(mp.mpc(t,2*mp.pi),p,sigma)
        error=abs(continued-reversed_value)
        limits.append(clean(dict(t=t,value=L,derivative=der,kms_continued=continued,reversed=reversed_value,
            kms_difference=error,passed=bool(error<mp.mpf('1e-50')),vacuum_commutator_norm=4*mp.im(L)**2)))
        for M in [4,16,64,256]:
            a=return_block(M,p,'.2',t,sigma,24);b=return_block(M,p,'.2',t,sigma,32);diff=abs(a['value']-b['value'])
            row=clean(dict(coarse=a,fine=b,difference=diff,rounding_allowance='1e-48',
                second_scaled=M*M*(b['value']-L-1j*der/M),passed=bool(diff<=a['tail']+b['tail']+mp.mpf('1e-48'))))
            rows.append(row)
            with raw.open('a',encoding='utf-8') as out:out.write(json.dumps(row)+'\n');out.flush()
            print('Gaussian',M,str(t),'vacuum',mp.nstr(b['vacuum_commutator_norm'],8),'passed',row['passed'],flush=True)
    for u in map(mp.mpf,['-2','-.4','.7','2']):
        err=abs(density(-u,p,sigma)/(mp.exp(2*mp.pi*u)*density(u,p,sigma))-1)
        checks.append(clean(dict(u=u,error=err,passed=bool(err<mp.mpf('1e-55')))))
    paths=['code/ds3_gaussian_field.py','code/check_ds3_gaussian_field.py','code/ds3_fock_vacuum.py','code/ds3_two_matter.py']
    data=dict(status='passed' if all(x['passed'] for x in limits+rows+checks) else 'failed',p=str(p),sigma=str(sigma),
        limits=limits,rows=rows,detailed_balance=checks,working_dps=60,
        multiplier_bound=clean(multiplier_bound(sigma)),block_bound=clean(multiplier_bound(sigma)*block_bound('.2')),
        scope='Real Gaussian local-field vacuum-return and limiting thermal spectrum; no new full U or fourpoint quadrature. KMS belongs to the heavy-observer limit, not finite mass.',
        dependencies=dict(python=platform.python_version(),mpmath=mp.__version__),
        hashes={q:hashlib.sha256((root/q).read_bytes()).hexdigest() for q in paths})
    dest.write_text(json.dumps(data,indent=2),encoding='utf-8');print(data['status'],flush=True)
    if data['status']!='passed':raise SystemExit(1)

if __name__=='__main__':main()
