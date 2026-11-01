"""Bind a bounded paired-kernel diagnostic, preserving earlier evidence."""
from pathlib import Path
import json,hashlib,zipfile
from datetime import datetime,timezone
from audit_paired_kernel import audit

def main():
    root=Path(__file__).resolve().parents[1];base=root/'results/stage4/paired-kernel'
    if (base/'evidence.json').exists():raise RuntimeError('Already bound')
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    priorpath=root/'results/stage4/weak-continuation/evidence.json'
    prior=json.loads(priorpath.read_text())
    for p,h in prior['hashes'].items():assert sha(root/p)==h,p
    known={r['name']:r['sha256'] for r in json.loads((root/'results/stage4/weak-continuation/preservation.json').read_text())['old_archives']}
    archives=[]
    for p in sorted((root/'snapshots').glob('*.zip')):
        h=sha(p);side=p.with_suffix('.sha256');assert side.exists() or p.name in known,p
        if side.exists():assert side.read_text().split()[0]==h,p
        if p.name in known:assert h==known[p.name],p
        archives.append(dict(name=p.name,sha256=h))
    assert len(archives)==56 and len(known)==55
    same=0;changed=[]
    allowed={'README.md','PLAN.md','NEXT.md','journal/2026-11-01.md','code/project_audit.py','results/project-audit.json','results/project-files.json'}
    with zipfile.ZipFile(root/'snapshots/2026-11-01-weak-continuation.zip') as z:
        for p in z.namelist():
            if p.endswith('/'):continue
            if (root/p).is_file() and sha(root/p)==hashlib.sha256(z.read(p)).hexdigest():same+=1
            else:assert p in allowed,p;changed.append(p)
    (base/'preservation.json').write_text(json.dumps(dict(status='passed',old_files_unchanged=same,
        changed_entries=changed,prior_bound_files=len(prior['hashes']),old_archives=archives),indent=2))
    sources={
        'sources/kim-scarinci-2112.13329.pdf':'https://arxiv.org/pdf/2112.13329v2',
        'sources/garoufalidis-kashaev-2306.01331v2.pdf':'https://arxiv.org/pdf/2306.01331v2',
        'sources/garoufalidis-kashaev-1706.08132v2.pdf':'https://arxiv.org/pdf/1706.08132v2',
        'sources/kashaev-beta-1403.1298v1.pdf':'https://arxiv.org/pdf/1403.1298v1'}
    (base/'source-record.json').write_text(json.dumps(dict(sources=[dict(file=p,url=u,sha256=sha(root/p)) for p,u in sources.items()],
        prior_evidence={priorpath.relative_to(root).as_posix():sha(priorpath)},
        reading_scope='Partial source reading: GK local-fields pages 4,5,8; GK index pages 30,31; Kashaev beta page 4 rendered and reviewed. KS scalar conventions and Proposition 6.15. Does not certify the paired pentagon.'),indent=2))
    replay=audit(root);assert replay['status']=='passed' and replay['checks']==1290,replay
    paths=['code/'+p for p in ['paired_kernel_heat.py','run_paired_kernel_heat.py',
           'check_paired_kernel.py','audit_paired_kernel.py','build_paired_kernel.ps1','record_paired_kernel.py']]
    paths+=['manuscript/paired-kernel/main.tex','output/pdf/de-sitter-paired-kernel.pdf',
            'notes/stage4-paired-kernel.md','literature/stage4-paired-kernel-sources.md']
    paths+=[p.relative_to(root).as_posix() for p in sorted(base.iterdir()) if p.is_file()]
    evidence=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),status='diagnostic stage bound; project audit follows',
        scientific_status='bounded-kernel formulation and empirical heat-tail study; pentagon not evaluated',pages=6,
        runs=240,nested_prefixes=True,diagnostic_checks=30,saved_arithmetic_replay_checks=1290,
        compilers=['pdfLaTeX twice exit 0','desktop success'],
        visual_review=dict(pages=list(range(1,7)),longest_edge_pixels=1100,status='all inspected'),
        hashes={p:sha(root/p) for p in paths})
    (base/'evidence.json').write_text(json.dumps(evidence,indent=2))
    print(json.dumps(dict(bound_files=len(paths),old_files_unchanged=same,old_archives=len(archives),scientific_status=evidence['scientific_status'])))
if __name__=='__main__':main()
