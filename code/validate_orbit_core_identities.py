"""Independent high-precision branch and ODE diagnostics; not a proof of FP completeness."""
from pathlib import Path
from datetime import datetime, timezone
import argparse, hashlib, json, platform
import mpmath as mp

ROOT=Path(__file__).resolve().parent.parent

def scalar(z):
    z=mp.mpc(z)
    return {"real":mp.nstr(z.real,65),"imag":mp.nstr(z.imag,65)}

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--label',required=True)
    args=parser.parse_args()
    if not args.label.replace('-','').isalnum(): parser.error('plain alphanumeric/hyphen label required')
    dest=ROOT/'results/stage3/orbit-core'/f'{args.label}.json'
    if dest.exists(): parser.error('use a new label; previous diagnostics are immutable')
    mp.mp.dps=70
    tol=mp.mpf('1e-55')
    rows=[]
    def check(kind,params,lhs,rhs):
        err=abs(lhs-rhs)/max(mp.mpf(1),abs(lhs),abs(rhs))
        rows.append({'kind':kind,**params,'lhs':scalar(lhs),'rhs':scalar(rhs),
                     'scaled_error':mp.nstr(err,20),'tolerance':str(tol),'passed':bool(err<tol)})
    # The original regularized hypergeometric expression and its two-branch
    # continuation are evaluated separately. A common smooth prefactor is 1.
    for n in (0,1,2):
        h=mp.mpf(n)+mp.mpf('.5')
        for ss in ('.4','.8','2.3'):
            s=mp.mpf(ss)
            for eps in (-1,1):
                for tt in ('4','8'):
                    tau=mp.mpf(tt); t=eps*tau; x=mp.exp(-tau)
                    z=(1-1j*mp.sinh(t))/2
                    lhs=mp.cosh(t)**h*mp.hyp2f1(h-1j*s,h+1j*s,n+1,z)/mp.gamma(n+1)
                    rhs=mp.mpc(0)
                    for sig in (-1,1):
                        q=h-1j*sig*s
                        amp=2**(2*q-h)*mp.exp(-1j*eps*mp.pi*q/2)*mp.gamma(2j*sig*s)/(mp.gamma(h+1j*sig*s)*mp.gamma(mp.mpf('.5')+1j*sig*s))
                        shape=(1+x*x)**h*(1+1j*eps*x)**(-2*q)*mp.hyp2f1(q,mp.mpf('.5')-1j*sig*s,1-2j*sig*s,4j*eps*x/(1+1j*eps*x)**2)
                        rhs+=amp*mp.exp(1j*sig*s*tau)*shape
                    check('two-end-mode-connection',{'abs_n':n,'s':ss,'t':str(t)},lhs,rhs)
    def jost(t,s):
        return mp.exp(1j*s*t)*mp.hyp2f1(mp.mpf('.5'),mp.mpf('.5'),1-1j*s,(1-mp.tanh(t))/2)
    for ss in ('.1','.6','2'):
        s=mp.mpf(ss)
        T=mp.gamma(mp.mpf('.5')-1j*s)**2/(mp.gamma(1-1j*s)*mp.gamma(-1j*s))
        for tt in ('-3','0','3'):
            t=mp.mpf(tt); F=lambda u:jost(u,s)
            check('jost-differential-equation',{'s':ss,'t':tt},-mp.diff(F,t,2)+mp.sech(t)**2*F(t)/4,s*s*F(t))
        fm=lambda t:jost(-t,s); fp=lambda t:jost(t,s)
        W=fm(0)*mp.diff(fp,0)-mp.diff(fm,0)*fp(0)
        check('transmission-wronskian',{'s':ss},2j*s/W,T)
        check('transmission-modulus',{'s':ss},abs(T)**2,mp.tanh(mp.pi*s)**2)
    inputs=['code/validate_orbit_core_identities.py','sources/kolchmeyer-liu-2411.08090v1.html',
            'sources/dlmf-15.8.html','sources/germain-pusateri-rousset-2018.pdf']
    failed=sum(not row['passed'] for row in rows)
    output={'status':'passed' if failed==0 else 'failed','created_at_utc':datetime.now(timezone.utc).isoformat(),
            'python':platform.python_version(),'mpmath':mp.__version__,'working_dps':mp.mp.dps,
            'checks':len(rows),'failed':failed,'max_scaled_error':str(max(mp.mpf(row['scaled_error']) for row in rows)),
            'scope':'Finite-point diagnostics of both time-end connection branches, the radial Jost ODE, Wronskian sign, and transmission modulus. Not interval-certified and not a substitute for estimates, density, FP closure, or local-field descent proofs.',
            'hashes':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs},'details':rows}
    dest.parent.mkdir(parents=True,exist_ok=True)
    with dest.open('x',encoding='utf-8') as f: json.dump(output,f,indent=2); f.write('\n')
    print(json.dumps({k:output[k] for k in ('status','checks','failed','max_scaled_error')}))
    print(dest)
    raise SystemExit(bool(failed))

if __name__=='__main__':main()
