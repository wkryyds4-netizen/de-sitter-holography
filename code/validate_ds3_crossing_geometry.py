"""Independent chart-area and analytic puncture-bound diagnostics."""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from ds3_crossing_dictionary import six_charts
from ds3_crossing_tail import puncture_tail_bound
from validate_ds3_crossing_dictionary import encode

ROOT=Path(__file__).resolve().parents[1]


def run():
    mp.mp.dps=45
    areas=[]; checks=[];tails=[]
    for n in [16,32]:
        xs,ws=mp.gauss_quadrature(n,'legendre');value=0
        for x,wx in zip(xs,ws):
            r=(x+1)/2;lo=mp.acos(r/2)
            for y,wy in zip(xs,ws):
                th=lo+(mp.pi-lo)*(y+1)/2;z=r*mp.exp(1j*th)
                value+=wx*wy*r*(mp.pi-lo)/8*sum(
                    jac/(1+abs(w)**2)**2 for w,jac in six_charts(z))
        err=abs(value/(mp.pi/4)-1)
        areas.append(dict(order=n,value=value,expected=mp.pi/4,error=err))
        checks.append(dict(name=f'area {n}',passed=bool(err<mp.mpf('1e-14'))))
    for z in [mp.mpc('.2','.3'),mp.mpc('-.4','.6'),mp.mpc('1.2','-.7'),mp.mpc('4','3')]:
        count=sum(abs(w)<1 and abs(1-w)>1 for w,_ in six_charts(z))
        checks.append(dict(name=f'unique chart {z}',passed=count==1))
    for ep in ['.001','.0001','1e-8','1e-12']:
        delta=mp.mpf(ep)/(1-mp.mpf(ep));c=128*mp.pi**3/3
        direct=c*mp.quad(lambda r:(2+r)*(12+mp.log(1/r))**2,[0,delta])
        bound=puncture_tail_bound(ep)
        err=abs(direct/bound-1)
        tails.append(dict(epsilon=ep,bound=bound,direct=direct,error=err))
        checks.append(dict(name=f'tail primitive {ep}',passed=bool(err<mp.mpf('1e-35'))))
    paths=['code/validate_ds3_crossing_geometry.py','code/ds3_crossing_tail.py',
           'code/ds3_crossing_dictionary.py']
    data=dict(dps=45,areas=areas,tails=tails,checks=checks,
              passed=sum(c['passed'] for c in checks),total=len(checks),
              hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths})
    (ROOT/'results/stage4/ds3-crossing-dictionary/geometry.json').write_text(
        json.dumps(encode(data),indent=2),encoding='utf-8')
    print(data['passed'],data['total']);print(encode(areas));print(encode(tails))
    if data['passed']!=data['total']:raise SystemExit(1)


if __name__=='__main__':run()
