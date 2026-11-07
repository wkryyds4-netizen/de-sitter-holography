"""Supplement: real 50-decimal check after discovering Windows longdouble=double.

Reads sealed first output; does not rerun or change that successful driver.
"""
from decimal import Decimal as D, localcontext
from pathlib import Path
from datetime import datetime, timezone
import hashlib,json,platform,sys

ROOT=Path(__file__).resolve().parents[1]
STAGE='susskind-thermal-chord-window'
SRC=f'results/{STAGE}/window-checks.json'
OUT=ROOT/f'results/{STAGE}/precision-checks.json'
if OUT.exists():raise FileExistsError(OUT)
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
assert sha(SRC)=='d36161666759fda1bfe351a0a09be636692931705fef7facdecf27a10a4cd6b0'
source=json.loads((ROOT/SRC).read_text())
rows=[]
with localcontext() as ctx:
 ctx.prec=50
 for row in source['rows']:
  m=D(int(row['m']));eps=1/m**3;cut=row['fine_cutoff'];order=row['taylor']['order']
  q=(-eps).exp();den=1-q;t=m.sqrt()
  qi=[(1-q**n)/den for n in range(cut+3)]
  off=[v.sqrt() for v in qi]
  term=[D(0)]*(cut+1);term[0]=(-m).exp();vec=term.copy()
  for k in range(1,order+1):
   nxt=[D(0)]*(cut+1);factor=t/k
   for n in range(min(k,cut)+1):
    v=(off[n]*term[n-1] if n else D(0))
    if n<cut:v+=off[n+1]*term[n+1]
    nxt[n]=factor*v
   term=nxt
   vec=[a+b for a,b in zip(vec,term)]
  norm2=sum(a*a for a in vec);norm=norm2.sqrt()
  psi=[a/norm*(-1 if n%2 else 1) for n,a in enumerate(vec)]
  probs=[v*v for v in psi];mean=sum(D(n)*v for n,v in enumerate(probs))
  chi=[(-m/2).exp()]
  for n in range(1,cut+1):chi.append(-t/D(n).sqrt()*chi[-1])
  dist=sum((a-b)**2 for a,b in zip(psi,chi)).sqrt()
  hmean=2*sum(off[n]*psi[n]*psi[n-1] for n in range(1,cut+1))
  r=q**2;omr=1-r;omr2=1-r*r;ps=psi+[D(0)]
  defect2=D(0);edeffect2=D(0);b2=D(0);bp=D(0)
  for n,v in enumerate(ps):
   c=(-eps*D(n+2)/2).exp();root=(1-c*c).sqrt()
   g=(qi[n]*omr2).sqrt()
   be=qi[n]*omr2/(2*c*root)*v;bf=-r*g/(2*c*root)*v
   prev=ps[n-1] if n else D(0)
   nxt=ps[n+1] if n+1<len(ps) else D(0)
   qe=omr/(2*eps.sqrt())*(off[n]*prev+off[n+1]*nxt)
   qf=-omr2.sqrt()/(2*eps.sqrt())*prev
   ee=((1-q*q)/den-omr*omr*qi[n]/2)/c*v
   ef=-r*g/(2*c)*v
   defect2+=(be+qe)**2+(bf+qf)**2
   edeffect2+=(ee-2*v)**2+ef**2
   b2+=be*be+bf*bf
   bp+=2*bf*(-g/(2*root)*v)
  values=dict(seed_overlap=psi[0],number_mean=mean,coherent_distance=dist,
      scaled_norm_squared=norm2,boost_sum_norm=defect2.sqrt(),
      E_eigenvector_defect=edeffect2.sqrt(),B_second=b2,BP_commutator_expectation_imag=bp)
  reference=dict(seed_overlap=row['seed_overlap'],number_mean=row['number_mean'],
      coherent_distance=row['coherent_distance'],scaled_norm_squared=row['taylor']['scaled_norm_squared'],
      boost_sum_norm=row['observables']['2']['boost_sum_norm'],
      E_eigenvector_defect=row['observables']['2']['E_eigenvector_defect'],
      B_second=row['observables']['2']['B']['second'],
      BP_commutator_expectation_imag=row['observables']['2']['BP_commutator_expectation_imag'])
  rel={k:float(abs(v-D(str(reference[k])))/abs(v)) for k,v in values.items()}
  assert max(rel.values())<2e-10
  rows.append(dict(m=int(m),decimal_precision=ctx.prec,cutoff=cut,order=order,
     values={k:str(v) for k,v in values.items()},relative_difference=rel,
     h_mean=str(hmean),h_mean_over_sqrt_m=str(hmean/t)))
  print(json.dumps(dict(m=int(m),max_relative=max(rel.values()))),flush=True)
result=dict(status='passed',created_at_utc=datetime.now(timezone.utc).isoformat(),rows=rows,
 correction='On this Windows host numpy.longdouble has nmant=52. The first output compares two double-precision evaluation paths, not 53/64-bit significands. Its limitations string overstates this precision check and is superseded here. No sealed bytes changed. Decimal uses 50 significant decimal digits.',
 scope='Positive Jacobi Taylor recurrence at actual 50-decimal precision, compared against sealed double-precision observables for d=2. Analytic remainder bounds are exact-arithmetic, no interval-rounding certification.',
 input_sha256={SRC:sha(SRC),'code/check_thermal_chord_precision.py':sha('code/check_thermal_chord_precision.py')},
 environment=dict(python=sys.version,platform=platform.platform(),decimal_precision=50))
with OUT.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(dict(status='passed',sha256=sha(OUT.relative_to(ROOT)),
 max_relative=max(v for r in rows for v in r['relative_difference'].values()))))
