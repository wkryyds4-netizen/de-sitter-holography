"""Freeze a diagnostic stage without converting its unresolved claim to success."""
from pathlib import Path
import hashlib,json,zipfile
from datetime import datetime,timezone
from audit_unitary_flip import audit
def main():
    root=Path(__file__).resolve().parents[1];base=root/'results/stage4/unitary-flip'
    if (base/'evidence.json').exists():raise RuntimeError('Already frozen')
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    prior=json.loads((root/'results/stage4/mass-filter/evidence.json').read_text())
    for p,h in prior['hashes'].items():assert sha(root/p)==h,p
    known={r['name']:r['sha256'] for r in json.loads((root/'results/stage4/mass-filter/preservation.json').read_text())['old_archives']}
    archives=[]
    for p in sorted((root/'snapshots').glob('*.zip')):
        h=sha(p);side=p.with_suffix('.sha256');assert side.exists() or p.name in known,p
        if side.exists():assert side.read_text().split()[0]==h,p
        if p.name in known:assert known[p.name]==h,p
        archives.append(dict(name=p.name,sha256=h))
    assert len(archives)==53 and len(known)==52
    same=0;changed=[]
    allowed={'README.md','PLAN.md','NEXT.md','code/project_audit.py','results/project-audit.json','results/project-files.json'}
    with zipfile.ZipFile(root/'snapshots/2026-10-31-mass-filter.zip') as z:
        for p in z.namelist():
            if p.endswith('/'):continue
            if (root/p).is_file() and sha(root/p)==hashlib.sha256(z.read(p)).hexdigest():same+=1
            else:assert p in allowed,p;changed.append(p)
    (base/'preservation.json').write_text(json.dumps(dict(status='passed',old_files_unchanged=same,
        changed_entries=changed,prior_bound_files=len(prior['hashes']),old_archives=archives),indent=2))
    sources={
        'sources/chandra-liska-post-2610.30247v1.html':'https://arxiv.org/html/2610.30247v1',
        'sources/kim-trilogy-2308.02934v2.html':'https://arxiv.org/html/2308.02934v2',
        'sources/kim-scarinci-2112.13329.pdf':'https://arxiv.org/pdf/2112.13329v2'}
    (base/'source-record.json').write_text(json.dumps(dict(sources=[dict(file=p,url=u,sha256=sha(root/p)) for p,u in sources.items()],
        online_only=['https://arxiv.org/html/2501.01486v1'],
        version_note='KS downloaded from unversioned pdf endpoint; title page confirms v2, 21 June 2024. Prop.6.15 in this PDF, (3.8) in Kim v2.',
        prior_evidence={'results/stage4/mass-filter/evidence.json':sha(root/'results/stage4/mass-filter/evidence.json')}),indent=2))
    a=audit(root);assert a['status']=='passed' and a['checks']==161,a
    code=['unitary_flip.py','run_unitary_flip.py','check_unitary_flip.py','audit_unitary_flip.py','record_unitary_flip.py','build_unitary_flip.ps1']
    paths=['code/'+p for p in code]+[p.relative_to(root).as_posix() for p in sorted((root/'code/unitary-flip-prototypes').glob('*.py'))]
    paths+=['manuscript/unitary-flip/main.tex','output/pdf/de-sitter-unitary-flip.pdf','notes/stage4-unitary-flip.md','literature/stage4-unitary-flip-sources.md']
    paths+=[p.relative_to(root).as_posix() for p in sorted(base.iterdir()) if p.is_file()]
    data=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),status='diagnostic stage bound; project audit follows',
        scientific_status='non-flat pentagon reproduction unresolved; no continuum counterexample',pages=5,
        matched_protocol_runs=21,scalar_algebra_checks=36,evidence_replay_checks=161,
        compilers=['pdfLaTeX twice exit 0','desktop success'],
        visual_review=dict(pages=[1,2,3,4,5],longest_edge_pixels=1100,status='all inspected'),
        hashes={p:sha(root/p) for p in paths})
    (base/'evidence.json').write_text(json.dumps(data,indent=2))
    print(json.dumps(dict(bound_files=len(paths),old_files_unchanged=same,old_archives=len(archives),scientific_status=data['scientific_status'])))
if __name__=='__main__':main()
