"""Independent diagnostics for stage3-clock-rate.md; not interval certification.

Old bound code is imported read-only. Time arguments are the two proper-time
pairing differences (u,v), not an implicit four-operator display order.
"""
from __future__ import annotations
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib, json, platform
from datetime import datetime, timezone
import numpy as np
import scipy
import mpmath as mp
from quantized_observer import DoublePacket, right_kernel, composite_gauss
from semiclassical_clock import stationary_transition

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'results'/'stage3'/'clock-rate'
CHECKS = []
def check(name, value, limit):
    value = float(value)
    CHECKS.append(dict(name=name, value=value, limit=float(limit), passed=bool(value <= limit)))
def complex_json(z):
    return [float(np.real(z)), float(np.imag(z))]
def gamma4(c, nu, x):
    return abs(c.gamma(c.mpf('.25')+.5j*(nu+x)))**2 * abs(c.gamma(c.mpf('.25')+.5j*(nu-x)))**2

def exact_ratio(nu, p, q, M, dps):
    c = mp.mp.clone(); c.dps = dps
    nu,p,q,M = map(c.mpf, map(str, (nu,p,q,M)))
    m = M+p/(2*c.pi); mm = M+q/(2*c.pi)
    a = c.sqrt(mm*mm-c.mpf('.25')); b = c.sqrt(m*m-c.mpf('.25'))
    delta = mm-m; d = a-b; S = a+b
    v = lambda x: c.exp(-c.pi*x)*gamma4(c,nu,x)/(16*c.pi**3)
    direct = mm*c.tanh(c.pi*a)*c.cosh(c.pi*b)**2*gamma4(c,nu,d)*gamma4(c,nu,S)/(16*c.pi**5*v(delta))
    T = S*c.exp(c.pi*S)*gamma4(c,nu,S)/(8*c.pi**2)
    factored = 2*mm/S*c.tanh(c.pi*a)*(1+c.exp(-2*c.pi*b))**2*T*v(d)/v(delta)
    return dict(direct=direct, factored=factored, delta=delta, local_m=m)

class FlatTask:
    def __init__(self, M, order=24, p_order=40, extra=14., nu=1., sigma=.5):
        packet = DoublePacket()
        p, wp = packet.nodes(p_order)
        self.outer = wp*packet.probability(p)
        self.x = p/(2*np.pi)
        self.y, wy = composite_gauss(0., packet.intervals[-1][1]/(2*np.pi)+nu+extra, order)
        self.de = self.y[None,:]-self.x[:,None]
        self.dg = self.y[None,:]-self.y[:,None]
        self.ext0 = stationary_transition(self.de,nu)*np.exp(-(sigma*self.de)**2)*wy
        self.grid0 = stationary_transition(self.dg,nu)*np.exp(-(sigma*self.dg)**2)*wy
        self.ext1 = self.ext0*self.de/2
        self.grid1 = self.grid0*self.dg/2
        if M is not None:
            mx,my = M+self.x,M+self.y
            sx,sy = np.sqrt(mx**2-.25),np.sqrt(my**2-.25)
            self.ext = right_kernel(nu,sy[None,:],sx[:,None])*(my/sy)[None,:]*np.exp(-(sigma*self.de)**2)*wy
            self.grid = right_kernel(nu,sy[None,:],sy[:,None])*(my/sy)[None,:]*np.exp(-(sigma*self.dg)**2)*wy
        else:
            self.ext,self.grid = self.ext0,self.grid0
    def values(self,u,v):
        E = self.ext*np.exp(1j*u*self.de)
        H_u,H_v = E.sum(axis=1),(self.ext*np.exp(1j*v*self.de)).sum(axis=1)
        G_v = (self.grid*np.exp(1j*v*self.dg)).sum(axis=1)
        E0 = self.ext0*np.exp(1j*u*self.de)
        E1 = self.ext1*np.exp(1j*u*self.de)
        H0u,H0v = E0.sum(axis=1),(self.ext0*np.exp(1j*v*self.de)).sum(axis=1)
        H1u,H1v = E1.sum(axis=1),(self.ext1*np.exp(1j*v*self.de)).sum(axis=1)
        G0v = (self.grid0*np.exp(1j*v*self.dg)).sum(axis=1)
        G1v = (self.grid1*np.exp(1j*v*self.dg)).sum(axis=1)
        vals = dict(C2=self.outer@H_u,AABB=self.outer@(H_u*H_v),ABBA=self.outer@(E@G_v))
        coefs = dict(C2=self.outer@H1u,AABB=self.outer@(H1u*H0v+H0u*H1v),ABBA=self.outer@(E1@G0v+E0@G1v))
        return vals,coefs

def main():
    point_rows=[]
    for nu,p,q in [(1,3,9),(1,9,3),(.4,0,2),(2.3,80,0),(1,3,3)]:
        values=[]
        for M in (40,80,160,320):
            row=exact_ratio(nu,p,q,M,60)
            check(f'exact-factorization:{nu,p,q,M}',abs(row['direct']-row['factored']),1e-48)
            target=row['delta']/2
            values.append(dict(M=M,ratio=str(row['direct']),coefficient_target=str(target),scaled_difference=str(M*(row['direct']-1)),coefficient_error=float(abs(M*(row['direct']-1)-target))))
        high=exact_ratio(nu,p,q,320,85)
        precision_ctx=mp.mp.clone(); precision_ctx.dps=85
        check(f'precision-refinement:{nu,p,q}',abs(precision_ctx.mpf(values[-1]['ratio'])-high['direct']),1e-48)
        # Selected fixed increments: first-coefficient error must decrease.
        check(f'coefficient-convergence:{nu,p,q}',values[-1]['coefficient_error']/max(values[-2]['coefficient_error'],1e-100),.7)
        point_rows.append(dict(nu=nu,p=p,q=q,values=values))
    # Large and M-dependent labels exercise exact global factorization, not a
    # claim that the fixed-p coefficient holds for drifting state families.
    global_rows=[]
    for p,q,M in [(0,10000,1),(10000,0,1),(1e6,1e6+3,20),(1e6,1,20)]:
        row=exact_ratio(1,p,q,M,70)
        check(f'global-factorization:{p,q,M}',abs(row['direct']-row['factored']),1e-52)
        global_rows.append(dict(p=p,q=q,M=M,ratio=str(row['direct'])))
    times=[(0.,0.),(.7,-.4),(1.,1.)]
    baseline=FlatTask(None)
    reference=[baseline.values(*uv) for uv in times]
    rows=[]
    for M in (10.,20.,40.,80.,160.):
        task=FlatTask(M)
        for uv,(ref,coef) in zip(times,reference):
            vals,_=task.values(*uv)
            for kind in ('C2','AABB','ABBA'):
                residual=vals[kind]-ref[kind]-coef[kind]/M
                rows.append(dict(M=M,u=uv[0],v=uv[1],channel=kind,value=complex_json(vals[kind]),limit=complex_json(ref[kind]),first_coefficient=complex_json(coef[kind]),scaled_difference=complex_json(M*(vals[kind]-ref[kind])),residual=complex_json(residual),scaled_second_residual=complex_json(M*M*residual)))
    fine=FlatTask(160.,order=32,p_order=56,extra=16.)
    for uv in times:
        vals,_=fine.values(*uv)
        for kind in ('C2','AABB','ABBA'):
            coarse=next(r for r in rows if r['M']==160 and (r['u'],r['v'])==uv and r['channel']==kind)
            old=complex(*coarse['value'])
            check(f'quadrature-and-tail-refinement:{uv,kind}',abs(vals[kind]-old),2e-10)
            r80=next(r for r in rows if r['M']==80 and (r['u'],r['v'])==uv and r['channel']==kind)
            e80=abs(complex(*r80['residual']))
            e160=abs(complex(*coarse['residual']))
            check(f'first-corrected-residual-contraction:{uv,kind}',e160/max(e80,1e-100),.4)
    step=1e-4
    for uv,(ref,coef) in zip(times,reference):
        u,v=uv
        pu=baseline.values(u+step,v)[0];mu=baseline.values(u-step,v)[0]
        pv=baseline.values(u,v+step)[0];mv=baseline.values(u,v-step)[0]
        for kind in ('AABB','ABBA'):
            derivative=(pu[kind]-mu[kind]+pv[kind]-mv[kind])/(4j*step)
            check(f'time-derivative-coefficient:{uv,kind}',abs(derivative-coef[kind]),2e-9)
    files=[Path(__file__),ROOT/'code'/'quantized_observer.py',ROOT/'code'/'semiclassical_clock.py',ROOT/'notes'/'stage3-clock-rate.md']
    metadata=dict(utc=datetime.now(timezone.utc).isoformat(),python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,mpmath=mp.__version__,model='dS2; G->0; nu=1 and sigma=.5 in packet quadrature; Q at every insertion',packet='DoublePacket: centers 3,9; halfwidth 1; equal lobe probabilities',hashes={str(x.relative_to(ROOT)):hashlib.sha256(x.read_bytes()).hexdigest() for x in files},limits='Diagnostics only; no interval certification, no verification of every heavy-tail P, no crossed channel.')
    study=dict(metadata=metadata,pointwise=point_rows,global_factorization=global_rows,fourpoint_and_twopoint=rows)
    validation=dict(metadata=metadata,checks=CHECKS,summary=dict(total=len(CHECKS),passed=sum(x['passed'] for x in CHECKS),failed=sum(not x['passed'] for x in CHECKS)))
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'study.json').write_text(json.dumps(study,ensure_ascii=False,indent=2),encoding='utf-8')
    (OUT/'validation.json').write_text(json.dumps(validation,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(validation['summary']))
    for item in CHECKS:
        if not item['passed']:print(json.dumps(item))
    for row in rows:
        if row['M']==160 and row['u']==row['v']==0:print(json.dumps(row))
    if validation['summary']['failed']:raise SystemExit(1)
if __name__=='__main__':main()