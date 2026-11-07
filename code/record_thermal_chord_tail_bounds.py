"""One-shot analytic tail audit, no heat-vector recomputation."""
from pathlib import Path
from decimal import Decimal as D,localcontext
from datetime import datetime,timezone
import hashlib,json,math
ROOT=Path(__file__).resolve().parents[1]
STAGE='susskind-thermal-chord-window'
SRC=f'results/{STAGE}/window-checks.json'
OUT=ROOT/f'results/{STAGE}/tail-bounds.json'
if OUT.exists():raise FileExistsError(OUT)
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
source=json.loads((ROOT/SRC).read_text())
rows=[]
with localcontext() as ctx:
 ctx.prec=60
 def log_chernoff(mu,k):
  assert k>mu
  return -D(mu)+D(k)*(1+(D(mu)/D(k)).ln())
 for row in source['rows']:
  m=D(int(row['m']));delta=2/m
  k=row['taylor']['order'];threshold=(k+2)//2
  logt=log_chernoff(2*m,threshold)/2
  cut=row['cutoff'];logc=((D(2)/(1-delta)).ln()+log_chernoff(2*m,cut+1))/2
  assert logt<D('-80') and logc<D('-45')
  rows.append(dict(m=int(m),taylor_order=k,poisson_threshold=threshold,
    scaled_vector_taylor_tail_bound=str(logt.exp()),
    scaled_vector_taylor_log_bound=str(logt),
    normalized_cutoff_vector_bound=str(logc.exp()),
    correction_needed=row['taylor']['scaled_vector_tail_bound']==0.))
result=dict(status='passed',created_at_utc=datetime.now(timezone.utc).isoformat(),rows=rows,
 correction='The first output records zero for the m=512 finite-matrix Taylor bound because scipy poisson.logsf underflows. Zero is not an error certificate. Supersede all those Taylor bounds by the state-specific positive-walk bound here; no sealed driver or vector rerun.',
 proof='For positive Tq<=T1 entrywise, the omitted Taylor vector has norm squared <= ||exp(t T1)e0||^2 - ||sum_{j=0}^K(t T1)^j/j! e0||^2. Missing pairs of powers require total even degree 2l>=K+1. Dividing by exp(2m) bounds the scaled-vector squared error by P(Poisson(2m)>=ceil((K+1)/2)). Evaluate with log Chernoff <= -mu+k(1+log(mu/k)), k>mu. The same crossing-level argument gives normalized Jacobi-cutoff error squared <= 2 P(Poisson(2m)>=cutoff+1)/(1-delta). These bound exact arithmetic only.',
 input_sha256={SRC:sha(SRC),'code/record_thermal_chord_tail_bounds.py':sha('code/record_thermal_chord_tail_bounds.py')})
with OUT.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(dict(status='passed',sha256=sha(OUT.relative_to(ROOT)),rows=rows)))
