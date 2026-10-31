"""Replay scalar raw arithmetic and normalization evidence without kernel code."""
import hashlib,json,math
from pathlib import Path

def audit(root):
    root=Path(root);out=root/'results/stage4/mass-filter';errors=[];count=0
    def ck(ok,msg):
        nonlocal count
        count+=1
        if not ok:errors.append(msg)
    def mean(v):return math.fsum(v)/len(v)
    def se(v):
        m=mean(v);return math.sqrt(math.fsum((x-m)**2 for x in v)/(len(v)*(len(v)-1)))
    data={n:json.loads((out/(n+'.json')).read_text()) for n in ['study','packet','diagnostics']}
    for name,d in data.items():
        for p,h in d['hashes'].items():ck(hashlib.sha256((root/p).read_bytes()).hexdigest()==h,'hash '+p)
    d=data['diagnostics'];ck(d['status']=='passed' and len(d['checks'])==23,'diagnostic completeness')
    for c in d['checks']:ck(c['passed'] and c['error']<c['tolerance'],'diagnostic '+c['name'])
    for name,file,n,keys in [('study','runs.jsonl',64,['correlated','independent','difference']),('packet','packet-runs.jsonl',16,['packet','conditional','difference'])]:
        raw=[json.loads(x) for x in (out/file).read_text().splitlines()]
        ck(raw==data[name]['rows'] and len(raw)==n,'raw equality '+name)
        for i,r in enumerate(raw):
            ck(all(math.isfinite(x) for key in keys for x in r[key]),'finite '+str((name,i)))
            if name=='study':
                for k in range(2):ck(abs(r['correlated'][k]-r['independent'][k]-r['difference'][k])<2e-14,'paired difference '+str((i,k)))
                ck(abs(r['correlated'][0]-r['raw_correlated'][0]-(1-r['baseline_ratio']))<2e-14,'baseline '+str(i))
        for s in data[name]['summaries']:
            rs=[r for r in raw if all(r[k]==s[k] for k in (['delta','u','power'] if name=='study' else ['power']))]
            ck(len(rs)==8,'seed count')
            for key in keys:
                for k in range(2):
                    v=[r[key][k] for r in rs]
                    ck(abs(mean(v)-s[key][k])<1e-14,'mean replay '+key)
                    ck(abs(se(v)-s[key+'_se'][k])<1e-14,'se replay '+key)
    packet=data['packet'];p=packet['parameters']
    ck(p['G']==1/(4*p['kappa']) and p['kappa']==4000 and p['M']==1000,'fixed theory parameters')
    ck(0<packet['normalizer']<.25 and packet['normalizer_tail_bound']<1e-30,'normalizer tail scope')
    # Normalizer truncation is bounded; its quadrature estimate is not a certificate.
    independent=json.loads((root/'results/stage4/principal-eikonal/refinement.json').read_text())
    comparisons=[]
    for s in data['study']['summaries']:
        if s['power']!=21 or s['delta'] not in [.4,.8]:continue
        old=next(r for r in independent['summaries'] if (r['sigma'],r['gamma'],r['u'])==(s['delta'],.25,s['u']))
        for k in range(2):
            difference=abs(s['independent'][k]-old['mean'][k]);sem=math.hypot(s['independent_se'][k],old['empirical_stderr'][k])
            ck(difference<6*sem+1e-6,'previous protocol control')
            comparisons.append(dict(delta=s['delta'],u=s['u'],component=k,difference=difference,combined_se=sem))
    return dict(status='passed' if not errors else 'failed',checks=count,errors=errors,previous_protocol_comparisons=comparisons,
                scope='Raw arithmetic, protocol metadata and empirical controls; not a physical trace proof or total numerical error certificate.')

if __name__=='__main__':
    root=Path(__file__).resolve().parents[1];p=root/'results/stage4/mass-filter/audit.json'
    if p.exists():raise RuntimeError('Do not overwrite')
    result=audit(root);result['hashes']={'code/audit_mass_filter.py':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    p.write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
