"""Schwarz pairing and complete pure-recoil graph ratios, not physical OTOCs."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,platform,sys
import mpmath as mp
import sympy as sp
from recoil_internal import parameters,barnes
from recoil_vertex import triples,logF,logE
from recoil_bessel_joint import parse_complex,gauss_segment
from recoil_bessel_limit import J

ROOT=Path(__file__).resolve().parents[1]

def vertex_ratio(M,nu,omega,u=0,v=0):
    ts,ms=triples(M,nu,omega,u,v);ref,_=triples(M,nu,0,0,0)
    return mp.fprod(ms)/M**4*mp.exp(mp.fsum(logE(t) for t in ts)-mp.fsum(logE(t) for t in ref))

def vertex_pair_from_gamma(M,nu,omega,u,v):
    ts,_=triples(M,nu,omega,u,v);ref,_=triples(M,nu,0,0,0)
    lhs=mp.fsum(logF(t) for t in ts)
    sharp=mp.fsum(mp.conj(logF([mp.conj(x) for x in t])) for t in ts)
    norm=mp.fsum(mp.re(logF(t)) for t in ref)
    return mp.exp((lhs+sharp)/2-norm),mp.exp(mp.re(lhs)-norm)

def F(y):return mp.sqrt(1+y*y)-y*mp.asinh(y)

def auxiliary_window(K,time,order,dps):
    with mp.workdps(dps):
        K,time=map(mp.mpf,(K,time));lo=mp.mpf('2.5');hi=mp.mpf(3)
        endrate=2*(-F(hi)-1+time*hi)
        norm=abs(J(K,0))**2
        nodes,weights=mp.gauss_quadrature(order,'legendre')
        fn=lambda y:mp.exp(2*K*time*y-K*endrate)*abs(J(K,K*y))**2/norm
        value,_=gauss_segment(fn,lo,hi,nodes,weights,step=mp.mpf('.1'))
        endpoint=1/(K*2*(mp.asinh(hi)+time)*mp.sqrt(1+hi*hi))
        return value,dict(K=str(K),time=str(time),y_lower=str(lo),y_upper=str(hi),order=order,dps=dps,
            measure='dy, not d_alpha',maximum_step='.1',scaled_integral=str(value),extracted_exponent=str(K*endrate),
            endpoint_scaled_prediction=str(endpoint),endpoint_ratio=str(value/endpoint),
            reduced_saddle_rate=str(2*(mp.cosh(time)-1)),window_endpoint_rate=str(endrate),
            scope='auxiliary paired Bessel kernel only; exact graph window has an asymptotic theorem, not new quadrature here')

def run():
    mp.mp.dps=90;checks=[];vertices=[];symmetric=[];transverse=[];bessel=[]
    def check(name,error,tol):
        checks.append(dict(name=name,error=str(abs(error)),tolerance=str(tol),passed=bool(abs(error)<=tol)))
    M,w,u,v,t=sp.symbols('M w u v t',real=True)
    ms=sp.Matrix([M+w/2+u,M-w/2+v,M+w/2-u,M-w/2-v])
    check('mass-jacobian',int(abs(ms.jacobian([M,w,u,v]).det())-4),0)
    theta=[-sp.pi/2,sp.I*t,sp.pi/2,sp.pi+sp.I*t]
    shifts=[theta[0]-(theta[3]-2*sp.pi)]+[theta[i]-theta[i-1] for i in range(1,4)]
    check('external-transform-phase',int(sp.simplify(sum(shifts[i]*ms[i] for i in range(4))-2*sp.pi*M+2*sp.I*t*w)),0)
    for j,omega in enumerate([mp.mpf('.7'),mp.mpc('.4','.9'),mp.mpc(0,'1.3')]):
        M=mp.mpf(20);nu=mp.mpf(5);u=mp.mpf('.6');v=mp.mpf('1.1')
        _,ms=triples(M,nu,omega,u,v);_,rev=triples(M,nu,-omega,u,v)
        aa,rr=parameters([ms[0],nu,ms[1],ms[2],nu,ms[3]])
        ar,rs=parameters([rev[0],nu,rev[1],rev[2],nu,rev[3]])
        for z in [mp.mpc('.23','-22.7'),mp.mpc('-1.4','-26.2')]:
            check(f'original-kernel-even-{j}-{z}',barnes(z,aa,rr)/barnes(z,ar,rs)-1,mp.mpf('1e-75'))
        pair,wrong=vertex_pair_from_gamma(M,nu,omega,u,v)
        exact=vertex_ratio(M,nu,omega,u,v)
        check(f'holomorphic-vertex-pair-{j}',pair/exact-1,mp.mpf('1e-75'))
        check(f'vertex-Schwarz-{j}',vertex_ratio(M,nu,mp.conj(omega),u,v)/mp.conj(exact)-1,mp.mpf('1e-75'))
        vertices.append(dict(M=str(M),nu=str(nu),omega=str(omega),u=str(u),v=str(v),
            continued_pair=str(exact),gamma_pair=str(pair),wrong_new_modulus_square=str(wrong),
            relative_difference_wrong=str(abs(wrong/exact-1))))
    for K,alpha in [(4,mp.mpf('.7')),(8,mp.mpf('7.5')),(8,mp.mpf('15.5')),(16,mp.mpf('31.2'))]:
        K=mp.mpf(K);ak=K*mp.besselk(alpha,K);ai=K*mp.besseli(alpha,K)
        exact=abs(J(K,alpha))**2
        kk=ak*ak;ii=mp.pi**2*ai*ai;cross=2*mp.pi*mp.sin(mp.pi*alpha)*ak*ai
        check(f'Bessel-pair-interference-{K}-{alpha}',(kk+ii+cross-exact)/(ak+mp.pi*ai)**2,mp.mpf('1e-75'))
        bessel.append(dict(K=str(K),alpha=str(alpha),K_square=str(kk),I_square=str(ii),interference=str(cross),exact_square=str(exact)))
    old=json.loads((ROOT/'results/stage4/recoil-bessel-joint/diagnostics.json').read_text())
    new=json.loads((ROOT/'results/stage4/recoil-transverse/diagnostics.json').read_text())
    growth=json.loads((ROOT/'results/stage4/recoil-transverse/consequences.json').read_text())
    for row in old['joint_rows']:
        K,nu,alpha=map(mp.mpf,(row['K'],row['nu'],row['alpha']));M=nu*nu/(2*K);y=alpha/K;rho=K**3/nu**2
        ref=next(r for r in old['joint_rows'] if mp.mpf(r['K'])==K and mp.mpf(r['alpha'])==0)
        val=parse_complex(row['value']);v0=parse_complex(ref['value']);E=mp.mpf(row['envelope'])
        vr=vertex_ratio(M,nu,1j*alpha)
        paired=vr*abs(val/v0)**2;lead=vr*abs(J(K,alpha)/J(K,0))**2
        geometric=mp.exp(2*K*(F(y)-1))/mp.sqrt(1+y*y)
        boundscale=abs(vr)*E**2/abs(v0)**2
        with mp.workdps(110):high=vertex_ratio(M,nu,1j*alpha)
        check(f'complete-vertex-precision-{K}-{y}',abs(vr/high-1),mp.mpf('1e-75'))
        symmetric.append(dict(K=str(K),nu=str(nu),M=str(M),alpha=str(alpha),y=str(y),rho=str(rho),
            vertex_ratio=str(vr),complete_graph_ratio=str(paired),Bessel_pair_prediction=str(lead),
            geometric_prediction=str(geometric),ratio_to_geometric=str(paired/geometric),
            log_magnitude_over_2K=str(mp.log(abs(paired))/(2*K)),rate_prediction=str(abs(F(y))-1),
            normalized_pair_error=str(abs(paired-lead)/boundscale),pair_error_over_rho=str(abs(paired-lead)/(boundscale*rho)),
            provenance='saved exact-multiplier finite branches; new paired graph uses analytic continued vertex, not new full quadrature'))
    for g in growth['growth_rows']:
        K,nu=map(mp.mpf,(g['K'],g['nu']));M=nu*nu/(2*K);alpha=2*K;u=mp.mpf('.3')*nu;v=mp.mpf('.2')*nu
        vr=vertex_ratio(M,nu,1j*alpha,u,v)/vertex_ratio(M,nu,1j*alpha)
        pair=vr*mp.mpf(g['raw_internal_modulus_ratio'])**2
        asym=mp.mpf(g['asymptotic_ratio'])**2
        transverse.append(dict(K=str(K),nu=str(nu),alpha=str(alpha),u=str(u),v=str(v),
            analytic_vertex_ratio=str(vr),complete_transverse_ratio=str(pair),asymptotic_modulus=str(asym),
            relative_modulus_difference=str(abs(abs(pair)/asym-1)),observed_exponent=str(mp.log(abs(pair))/K),
            exponent_prediction=str(2*mp.mpf(g['delta']))))
    zero=old['zero_row'];z=parse_complex(zero['value']);ez=mp.mpf(zero['envelope'])
    zero_record=dict(K=zero['K'],alpha=zero['alpha'],nu=zero['nu'],paired_internal_over_E2=str(abs(z/ez)**2),
        scope='near an auxiliary Bessel cancellation; no division by J and no physical interpretation')
    yc=mp.findroot(F,(mp.mpf('1.4'),mp.mpf('1.6')));tc=mp.asinh(yc)
    check('crossover-equation',F(yc),mp.mpf('1e-80'))
    # A finite window of the auxiliary paired kernel already grows past the reduced saddle.
    windows=[]
    for ti in [0,1,mp.mpf('1.5')]:
        for yi in [2,3,4]:
            y=mp.mpf(yi);tt=mp.mpf(ti)
            windows.append(dict(t=str(tt),y=str(y),paired_exponent_per_K=str(2*(abs(F(y))-1+y*tt)),
                reduced_exponent_per_K=str(2*(F(y)-1+y*tt)),
                reduced_saddle_exponent_per_K=str(2*(mp.cosh(tt)-1))))
    window_integrals=[]
    for tt in [0,1]:
        for kk in [8,16,32]:
            low,_=auxiliary_window(kk,tt,20,55)
            high,row=auxiliary_window(kk,tt,32,75)
            err=abs(low/high-1)
            check(f'auxiliary-window-refinement-{kk}-{tt}',err,mp.mpf('1e-20'))
            row['relative_refinement']=str(err);window_integrals.append(row)
    inputs=['code/recoil_pairing.py','code/recoil_internal.py','code/recoil_vertex.py','code/recoil_bessel_limit.py',
        'code/recoil_bessel_joint.py','results/stage4/recoil-bessel-joint/diagnostics.json',
        'results/stage4/recoil-transverse/diagnostics.json','results/stage4/recoil-transverse/consequences.json',
        'sources/chandra-liska-post-2610.30247v1.html','sources/chen-stanford-tang-yang-2607.14042v1.html']
    data=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),status='passed' if all(c['passed'] for c in checks) else 'failed',
        scope='analytic chiral pairing, pure-recoil complete topology graph and limits of a single reduced exponent; not a physical OTOC',
        environment=dict(python=sys.version,platform=platform.platform(),mpmath=mp.__version__,sympy=sp.__version__,dps=[90,110]),
        checks=checks,vertex_rows=vertices,bessel_rows=bessel,symmetric_rows=symmetric,transverse_rows=transverse,
        auxiliary_zero=zero_record,crossover=dict(y=str(yc),time=str(tc)),window_rows=windows,window_integrals=window_integrals,
        hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs})
    out=ROOT/'results/stage4/recoil-pairing';out.mkdir(parents=True,exist_ok=True)
    (out/'diagnostics.json').write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=data['status'],checks=len(checks),failed=[c for c in checks if not c['passed']],crossover=data['crossover'],vertices=vertices,transverse=transverse),indent=2))
    if data['status']!='passed':raise SystemExit(1)
if __name__=='__main__':run()
