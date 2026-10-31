"""Check transverse multiplier scaling and compare saved complete branch evaluations."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib, json, platform, sys
import mpmath as mp
from recoil_transverse import effective,log_r,log_c
from recoil_bessel_joint import parse_complex
ROOT=Path(__file__).resolve().parents[1]

def run():
    mp.mp.dps=85; rows=[];checks=[];comparison=[]
    for uv in [('.3','.2'),('.4','-.3')]:
        for z in [mp.mpc('-1','.5'),mp.mpc('2','-.7')]:
            for ki in [16,64,256,1024]:
                K=mp.mpf(ki);nu=8*K**mp.mpf('1.75');u=nu*mp.mpf(uv[0]);v=nu*mp.mpf(uv[1]);t=K*z
                M,ws,k,lam=effective(K,nu,u,v);rho=K**3/nu**2
                residual=log_r(t,K,nu,u,v)-4j*lam*(t*t-t)
                with mp.workdps(110):
                    high=log_r(t,K,nu,u,v)-4j*effective(K,nu,u,v)[3]*(t*t-t)
                err=abs(residual-high)
                checks.append(dict(name=f'multiplier-precision-{uv}-{z}-{ki}',error=str(err),tolerance='1e-70',passed=bool(err<mp.mpf('1e-70'))))
                rows.append(dict(K=str(K),nu=str(nu),u_over_nu=uv[0],v_over_nu=uv[1],t_over_K=str(z),
                    k=str(k),lam=str(lam),rho=str(rho),phase_scale=str(K*K/nu),log_residual=str(residual),residual_over_rho=str(abs(residual)/rho)))
    new=json.loads((ROOT/'results/stage4/recoil-transverse/diagnostics.json').read_text())
    old=json.loads((ROOT/'results/stage4/recoil-bessel-joint/diagnostics.json').read_text())
    for r in new['joint_rows']:
        if r['family']!='fractional' or mp.mpf(r['alpha'])/mp.mpf(r['K'])!=2:continue
        K,nu,alpha,u,v=map(mp.mpf,(r['K'],r['nu'],r['alpha'],r['u'],r['v']))
        o=next(s for s in old['joint_rows'] if mp.mpf(s['K'])==K and mp.mpf(s['alpha'])==alpha)
        _,_,k,lam=effective(K,nu,u,v);q=k/K;y=alpha/K
        phi=lambda z:mp.sqrt(z*z+y*y)-y*mp.asinh(y/z)
        delta=phi(1)-phi(q);pref=((1+y*y)/(q*q+y*y))**mp.mpf('.25')
        cr=abs(mp.exp(log_c(K,nu,u,v)-log_c(K,nu,0,0)))
        actual=cr*abs(parse_complex(r['value'])/parse_complex(o['value']))
        asym=pref*mp.exp(K*delta)
        comparison.append(dict(K=str(K),nu=str(nu),y=str(y),u_over_nu=str(u/nu),v_over_nu=str(v/nu),q=str(q),
            delta=str(delta),C_modulus_ratio=str(cr),raw_internal_modulus_ratio=str(actual),
            asymptotic_ratio=str(asym),relative_difference=str(abs(actual/asym-1)),observed_exponent=str(mp.log(actual)/K),
            scope='saved finite branch evaluations plus exact C; analytic theorem controls omitted connectors asymptotically'))
    inputs=['code/transverse_consequences.py','code/recoil_transverse.py','code/recoil_bessel_joint.py',
        'results/stage4/recoil-transverse/diagnostics.json','results/stage4/recoil-bessel-joint/diagnostics.json']
    data=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),status='passed' if all(c['passed'] for c in checks) else 'failed',
        environment=dict(python=sys.version,platform=platform.platform(),mpmath=mp.__version__,dps=[85,110]),
        scope='multiplier precision checks and raw-growth samples; no numerical asymptotic threshold certificate',
        checks=checks,multiplier_rows=rows,growth_rows=comparison,
        hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs})
    (ROOT/'results/stage4/recoil-transverse/consequences.json').write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=data['status'],checks=len(checks),growth=comparison),indent=2))
    if data['status']!='passed':raise SystemExit(1)
if __name__=='__main__':run()
