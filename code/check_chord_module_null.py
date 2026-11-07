"""One-shot targeted null-generator covariance and positivity-form check."""
from pathlib import Path
from datetime import datetime,timezone
import json,hashlib,math
import numpy as np
import sympy as sp
from chord_module_frame import coherent_descendant,act,word
ROOT=Path(__file__).resolve().parents[1];STAGE='susskind-chord-module-frame'
SRC=f'results/{STAGE}/frame-checks.json';OUT=ROOT/f'results/{STAGE}/null-checks.json'
if OUT.exists():raise FileExistsError(OUT)
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
assert sha(SRC)=='f2227eca8b2c57f0d9292773bbf04046fb0c4687b8971eb98ab9625bc5615e27'
source=json.loads((ROOT/SRC).read_text());rows=[]
for rr in source['rows']:
 if (rr['j'],rr['k']) not in [(0,0),(2,3)]:continue
 m=rr['m'];eps=rr['epsilon'];M=rr['fine_cutoff']
 v=np.zeros((M+9,9),complex);v[:M+1,rr['k']]=coherent_descendant(m,rr['j'],M)
 tq=word(eps,2,'TQ',v);qt=word(eps,2,'QT',v);qv=act(eps,2,'Q',v)
 defect=np.linalg.norm(tq-qt-1j*qv)
 target=word(eps,2,'TQ',v,True)-word(eps,2,'QT',v,True)-1j*act(eps,2,'Q',v,True)
 assert np.linalg.norm(target)<5e-13
 rows.append(dict(m=m,j=rr['j'],k=rr['k'],null_commutator_defect=float(defect),
   relative=float(defect/np.linalg.norm(qv)),
   cancellation=float((np.linalg.norm(tq)+np.linalg.norm(qt)+np.linalg.norm(qv))/defect)))
exact=[]
for d in (1,2,3):
 M=7
 E=sp.diag(*[d+k for k in range(M+1)])
 L=sp.zeros(M+1)
 for k in range(M):L[k+1,k]=sp.sqrt((k+1)*(2*d+k))
 P=sp.I*(L-L.T)/2;Q=E-P
 A=sp.zeros(M+1)
 for k in range(M+1):
  A[k,k]=sp.sqrt(2*d+k)/sp.sqrt(2)
  if k<M:A[k,k+1]=sp.I*sp.sqrt(k+1)/sp.sqrt(2)
 assert sp.simplify(A.conjugate().T*A-Q)==sp.zeros(M+1)
 exact.append(dict(d=d,dimension=M+1,identity='Q=A^dagger A'))
result=dict(status='passed',created_at_utc=datetime.now(timezone.utc).isoformat(),
 rows=rows,exact_positivity_forms=exact,
 input_sha256={SRC:sha(SRC),'code/chord_module_frame.py':sha('code/chord_module_frame.py'),
  'code/check_chord_module_null.py':sha('code/check_chord_module_null.py')},
 scope='Finite moving vectors: [Btr,Qepsilon]-iQepsilon. Exact positivity form belongs to limiting relative Q, not finite-q LS matrices.',
 limitations='No finite-q positivity assertion, no opposite null generator or mass selection; main frame-checks contains environment and independent 50-decimal diagnostics.')
with OUT.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(dict(status='passed',sha256=sha(OUT.relative_to(ROOT)),rows=rows)))
