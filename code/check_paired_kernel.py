"""Check the source dictionary and heat controls; no pentagon is certified."""
from pathlib import Path
import json,hashlib,platform
import mpmath as m
import sympy as s

def main():
    root=Path(__file__).resolve().parents[1];base=root/'results/stage4/paired-kernel'
    out=base/'diagnostics.json'
    if out.exists():raise RuntimeError('Do not overwrite')
    checks=[];scalar=[];cylinder=[]
    def ck(name,ok,**kw):checks.append(dict(name=name,passed=bool(ok),**kw))
    h=s.symbols('h',positive=True);x,y,t=s.symbols('x y t',real=True)
    dot=x/(2*s.pi*s.sqrt(h));ddot=-s.sqrt(h)*y/s.pi
    ck('Gaussian-group quadratic form',s.simplify(2*s.pi*dot*ddot+x*y/s.pi)==0)
    ck('self-dual measure determinant',s.simplify(s.det(s.Matrix([dot,ddot]).jacobian([x,y]))+1/(2*s.pi**2))==0)
    a=s.symbols('a',positive=True)
    ck('vacuum heat norm',s.factor(1-2/(1+a)+1/(1+2*a)-2*a*a/((1+a)*(1+2*a)))==0)
    ck('chirp Gaussian average',s.simplify(1/(1+a)/(1+a/(1+a))-1/(1+2*a))==0)
    m.mp.dps=60
    def C(q,theta,index):
        z=m.exp(1j*theta);v=q**(1-index)
        return m.qp(-v*z,q*q)/m.qp(-v/z,q*q)
    encode=lambda z:[m.nstr(z.real,50),m.nstr(z.imag,50)]
    for hv in ['.5','1','2']:
        hh=m.mpf(hv);q=m.exp(-m.pi*hh);qt=m.exp(-m.pi/hh)
        for xx,yy in [('.2','.3'),('2','-4'),('3.2','-3.1')]:
            xx=m.mpf(xx);yy=m.mpf(yy);z=xx+1j*hh*yy
            phi=m.qp(-qt*m.exp(z/(1j*hh)),qt*qt)/m.qp(-q*m.exp(z),q*q)
            direct=phi/m.conj(phi)
            factored=C(qt,-xx/hh,hh*yy/m.pi)/C(q,hh*yy,xx/(m.pi*hh))
            inv=factored*C(qt,xx/hh,-hh*yy/m.pi)/C(q,-hh*yy,-xx/(m.pi*hh))
            target=m.exp(-1j*xx*yy/m.pi)
            row=dict(h=hv,x=str(xx),y=str(yy),direct=encode(direct),factored=encode(factored),
                inversion=encode(inv),target=encode(target),factor_error=m.nstr(abs(direct-factored),40),
                inversion_error=m.nstr(abs(inv-target),40));scalar.append(row)
            ck('continuous-index scalar factorization',abs(direct-factored)<m.mpf('1e-50'),h=hv,x=str(xx),y=str(yy))
            ck('paired inversion restoration',abs(inv-target)<m.mpf('1e-50'),h=hv,x=str(xx),y=str(yy))
    q=m.exp(-m.pi);theta=m.mpf('.7')
    for index in [-2,0,1,3]:
        val=C(q,theta,index)*C(q,-theta,-index);expected=m.exp(1j*index*theta)
        ck('integer cylinder inversion',abs(val-expected)<m.mpf('1e-50'),index=index)
    for label in ['.3','.5','1.2']:
        index=m.mpf(label);val=C(q,theta,index)*C(q,-theta,-index);target=m.exp(1j*index*theta)
        cylinder.append(dict(index=label,value=encode(val),branch_target=encode(target),
            difference=m.nstr(abs(val-target),40),target_monodromy=encode(m.exp(2j*m.pi*index))))
        ck('noninteger character is not single-valued',abs(m.exp(2j*m.pi*index)-1)>.1,index=label)
    rows=json.loads((base/'heat-study.json').read_text())['rows']
    control=[r for r in rows if r['h'] is None and r['power']==16]
    worst=max(abs(r['estimated_one_minus_C']/r['exact_one_minus_C']-1) for r in control)
    ck('analytic chirp control, empirical QMC tolerance',worst<.001,maximum_relative_deviation=worst)
    data=dict(status='passed' if all(c['passed'] for c in checks) else 'failed',checks=checks,
        scalar_rows=scalar,noninteger_rows=cylinder,precision=60,
        versions=dict(python=platform.python_version(),mpmath=m.__version__,sympy=s.__version__),
        scope='Dictionary and analytic Gaussian controls only; no certified QMC or pentagon evaluation.',
        hashes={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in
            ['code/check_paired_kernel.py','results/stage4/paired-kernel/heat-study.json']})
    out.write_text(json.dumps(data,indent=2));print(data['status'],len(checks))
if __name__=='__main__':main()
