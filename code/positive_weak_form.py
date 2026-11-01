"""The same moving-pole effect for one fixed normalized compact L2 state.

Its finite spectral first moment is a quadratic form, not <psi,A psi>
with A psi in L2. This is a mathematical test state, not a dS observer.
"""
from pathlib import Path
import json,hashlib
import mpmath as m
from weak_continuation import HC,point,nodes,normalized_direct
root=Path(__file__).resolve().parents[1];base=root/'results/stage4/weak-continuation'
out=base/'positive-form.json'
if out.exists():raise RuntimeError('Do not overwrite')
m.mp.dps=40;R=m.mpf('.04');power=8;mass=m.pi*R*R/9
def enc(z):return [m.nstr(z.real,40),m.nstr(z.imag,40)]
def form(H):
    offset=H-HC;total=[]
    for k in range(64):
        e=m.exp(2j*m.pi*k/64);p=(offset*m.conj(e)).real
        edge=-p+m.sqrt(p*p+R*R-abs(offset)**2)
        for s,weight in nodes(16,m.mp.dps):
            r=edge*s;u=H+r*e;w=point(u)[2]
            g=max(m.mpf(0),1-abs(u-HC)**2/R**2)**power
            total.append(weight*r*edge*g*normalized_direct(H,u)/(w*mass))
    return m.fsum(total)*2*m.pi/64
norm=m.fsum(weight*(1-s)**power for s,weight in nodes(16,m.mp.dps))*m.pi*R*R/mass
expected=m.pi/(point(HC)[2]*mass);rows=[]
for step in ['.002','.001','.0005']:
    d=m.mpf(step);values=[form(HC+d),form(HC-d),form(HC+1j*d),form(HC-1j*d)]
    deriv=(values[0]-values[1]+1j*(values[2]-values[3]))/(4*d)
    rows.append(dict(delta=step,values=[enc(z) for z in values],derivative=enc(deriv),
        expected=enc(expected),relative_error=m.nstr(abs(deriv/expected-1),30)))
data=dict(status='completed; finite-difference and quadrature errors empirical',precision=40,
    center=enc(HC),radius=str(R),power=8,mass=m.nstr(mass,40),norm=m.nstr(norm,40),rows=rows,
    state='psi(p(u))=sqrt(g(u)/(J(u)*mass)), positive real; g=(1-|u-Hc|^2/R^2)_+^8; normalized compact C3 state',
    normalized_form='integral |psi|^2 Phi^H(x+H y) dxdy / residue(H)',
    scope='Finite spectral first moment, not action of the unbounded factor on this vector; no physical dS-state claim.',
    hashes={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in ['code/positive_weak_form.py','code/weak_continuation.py']})
out.write_text(json.dumps(data,indent=2));print('complete',len(rows),data['norm'])
