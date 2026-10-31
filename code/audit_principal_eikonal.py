"""Replay saved raw estimates with Python arithmetic, without the evaluator."""
import hashlib,json,math
from pathlib import Path

def audit(root):
    root=Path(root);out=root/'results/stage4/principal-eikonal';errors=[];count=0
    def require(ok,msg):
        nonlocal count
        count+=1
        if not ok:errors.append(msg)
    def mean(x):return math.fsum(x)/len(x)
    def se(x):
        m=mean(x);return math.sqrt(math.fsum((v-m)**2 for v in x)/(len(x)*(len(x)-1)))
    data={name:json.loads((out/(name+'.json')).read_text()) for name in ['study','diagnostics','refinement','late']}
    for name,record in data.items():
        for p,h in record['hashes'].items():require(hashlib.sha256((root/p).read_bytes()).hexdigest()==h,'hash '+p)
    for name,n in [('diagnostics',60),('refinement',14)]:
        d=data[name];require(d['status']=='passed' and len(d['checks'])==n,'check count '+name)
        for c in d['checks']:require(c['passed'] and math.isfinite(c['error']) and c['error']<c['tolerance'],'criterion '+c['name'])
    for name,file,n in [('study','runs.jsonl',64),('refinement','refinement-runs.jsonl',32),('late','late-runs.jsonl',48)]:
        raw=[json.loads(line) for line in (out/file).read_text().splitlines()]
        require(raw==data[name]['rows'] and len(raw)==n,'raw consistency '+name)
        for i,r in enumerate(raw):
            require(all(math.isfinite(v) for key in ['value','raw_value','absolute_correction'] for v in r[key]),'finite '+str((name,i)))
            require(abs((r['value'][0]-r['raw_value'][0])-(1-r['baseline_ratio']))<2e-14,'baseline control '+str((name,i)))
            require(r['value'][1]==r['raw_value'][1],'imaginary control '+str((name,i)))
            require(r['absolute_raw']+1e-14>=math.hypot(*r['raw_value']),'absolute mass '+str((name,i)))
    for name,field,semfield in [('refinement','mean','empirical_stderr'),('late','raw_mean','raw_empirical_stderr')]:
        for s in data[name]['summaries']:
            rows=[r for r in data[name]['rows'] if all(r[k]==s[k] for k in (['sigma','gamma','u'] if name=='refinement' else ['u','rate','power']))]
            require(len(rows)==8,'scrambles '+str(s))
            for k in range(2):
                values=[r['value' if name=='refinement' else 'raw_value'][k] for r in rows]
                require(abs(mean(values)-s[field][k])<1e-14,'mean replay '+name)
                require(abs(se(values)-s[semfield][k])<1e-14,'stderr replay '+name)
    late=data['late']['summaries'];comparisons=[]
    for u in [4.,6.,8.]:
        a,b=[s for s in late if s['u']==u]
        for k in range(2):
            delta=abs(a['raw_mean'][k]-b['raw_mean'][k]);sem=math.hypot(a['raw_empirical_stderr'][k],b['raw_empirical_stderr'][k])
            require(delta<6*sem+1e-6,'late empirical consistency '+str((u,k)))
            comparisons.append(dict(u=u,component=k,difference=delta,combined_stderr=sem))
    require(late[1]['raw_mean'][1]>6*late[1]['raw_empirical_stderr'][1],'u4 resolved empirical positive sign')
    require(abs(late[-1]['raw_mean'][1])<2*late[-1]['raw_empirical_stderr'][1],'u8 remains unresolved')
    return dict(status='passed' if not errors else 'failed',checks=count,errors=errors,late_comparisons=comparisons,
                scope='Saved arithmetic, metadata and empirical consistency; not quadrature certification or new gravity evidence.')

if __name__=='__main__':
    root=Path(__file__).resolve().parents[1];result=audit(root)
    path=root/'results/stage4/principal-eikonal/audit.json'
    if path.exists():raise RuntimeError('Do not overwrite audit')
    result['hashes']={'code/audit_principal_eikonal.py':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    path.write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
