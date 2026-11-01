"""Save literal weak-integral checks without overwriting an earlier run."""
from pathlib import Path
import json,time,hashlib,platform
import mpmath as m
from weak_continuation import HC,constants,weak_subtracted,pole_centered_weak,fixed_integral

root=Path(__file__).resolve().parents[1];base=root/'results/stage4/weak-continuation'
base.mkdir(parents=True,exist_ok=True);raw=base/'runs.jsonl'
if raw.exists() or (base/'study.json').exists():raise RuntimeError('Do not overwrite')
m.mp.dps=40;radius=m.mpf('.04');power=8;rows=[]
def enc(z):return [m.nstr(z.real,40),m.nstr(z.imag,40)]
def emit(row):
    rows.append(row)
    with raw.open('a',encoding='utf-8') as f:f.write(json.dumps(row)+'\n')
    print(row['kind'],row.get('offset',row.get('delta',row.get('loop_nodes'))),flush=True)

for offset in [('0','0'),('.01','0'),('0','.01'),('.015','.01'),('-.01','.015'),('.06','0')]:
    H=HC+m.mpc(*offset);t=time.time()
    coarse=weak_subtracted(H,8,32,radius,power)
    fine=weak_subtracted(H,16,64,radius,power)
    direct=(pole_centered_weak(H,16,64,radius,power) if abs(H-HC)<radius
            else fixed_integral(H,16,64,radius,power,False))
    emit(dict(kind='integral',offset=offset,H=enc(H),subtracted_coarse=enc(coarse),subtracted_fine=enc(fine),direct=enc(direct),
              quadrature_difference=m.nstr(abs(fine-coarse),30),direct_difference=m.nstr(abs(fine-direct),30),
              direct_method='pole-centered' if abs(H-HC)<radius else 'fixed-disk',seconds=time.time()-t))

for delta in ['.004','.002','.001','.0005']:
    d=m.mpf(delta);points=[HC+d,HC-d,HC+1j*d,HC-1j*d]
    values=[weak_subtracted(H,8,32,radius,power) for H in points]
    result=(values[0]-values[1]+1j*(values[2]-values[3]))/(4*d)
    row=dict(kind='wirtinger',delta=delta,values=[enc(z) for z in values],derivative=enc(result),
             expected=enc(m.pi+0j),relative_error=m.nstr(abs(result/m.pi-1),30))
    if delta in ['.001','.0005']:
        independent=[pole_centered_weak(H,16,64,radius,power) for H in points]
        replay=(independent[0]-independent[1]+1j*(independent[2]-independent[3]))/(4*d)
        row.update(direct_values=[enc(z) for z in independent],direct_derivative=enc(replay),
                   direct_difference=m.nstr(abs(result-replay),30))
    emit(row)

loop_radius=m.mpf('.06');mass=m.pi*radius**2/(power+1);expected=2j*m.pi*mass
for loop_nodes,nr,nt in [(8,12,32),(16,16,48),(32,16,64)]:
    t=time.time();samples=[];total=m.mpc(0)
    for k in range(loop_nodes):
        e=m.exp(2j*m.pi*k/loop_nodes);H=HC+loop_radius*e
        value=fixed_integral(H,nr,nt,radius,power,False)
        total+=value*1j*loop_radius*e*2*m.pi/loop_nodes;samples.append(enc(value))
    emit(dict(kind='morera',loop_nodes=loop_nodes,source_radial=nr,source_angular=nt,loop_radius=str(loop_radius),
              samples=samples,integral=enc(total),expected=enc(expected),relative_error=m.nstr(abs(total/expected-1),30),seconds=time.time()-t))

data=dict(status='completed; quadrature empirical',precision=40,center=enc(HC),radius=str(radius),power=power,
          test='f(p(u))=(y(u)+pi*i)*g(u)/J(u); g=(1-|u-Hc|^2/radius^2)_+^8; compact C7 complex test, not a physical probability density',
          normalized_pairing='Jweak(H)=integral f(x,y)*Phi^H(x+H*y) dxdy / residue(H)',
          scope='Literal single-chiral weak pairing; not a bounded phase identity or finite-G state.',
          versions=dict(python=platform.python_version(),mpmath=m.__version__),rows=rows,
          hashes={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in ['code/weak_continuation.py','code/run_weak_continuation.py']})
(base/'study.json').write_text(json.dumps(data,indent=2));print('complete',len(rows),flush=True)
