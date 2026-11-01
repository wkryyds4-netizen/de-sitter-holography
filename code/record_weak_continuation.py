"""Bind a literal weak-continuation diagnostic, preserving earlier evidence."""
from pathlib import Path
import json,hashlib,zipfile
from datetime import datetime,timezone
from audit_weak_continuation import audit

def main():
    root=Path(__file__).resolve().parents[1];base=root/'results/stage4/weak-continuation'
    if (base/'evidence.json').exists():raise RuntimeError('Already bound')
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    priorpath=root/'results/stage4/flip-domain/evidence.json'
    prior=json.loads(priorpath.read_text())
    for p,h in prior['hashes'].items():assert sha(root/p)==h,p
    known={r['name']:r['sha256'] for r in json.loads((root/'results/stage4/flip-domain/preservation.json').read_text())['old_archives']}
    archives=[]
    for p in sorted((root/'snapshots').glob('*.zip')):
        h=sha(p);side=p.with_suffix('.sha256');assert side.exists() or p.name in known,p
        if side.exists():assert side.read_text().split()[0]==h,p
        if p.name in known:assert h==known[p.name],p
        archives.append(dict(name=p.name,sha256=h))
    assert len(archives)==55 and len(known)==54
    same=0;changed=[]
    allowed={'README.md','PLAN.md','NEXT.md','journal/2026-11-01.md','code/project_audit.py','results/project-audit.json','results/project-files.json'}
    with zipfile.ZipFile(root/'snapshots/2026-11-01-flip-domain.zip') as z:
        for p in z.namelist():
            if p.endswith('/'):continue
            if (root/p).is_file() and sha(root/p)==hashlib.sha256(z.read(p)).hexdigest():same+=1
            else:assert p in allowed,p;changed.append(p)
    (base/'preservation.json').write_text(json.dumps(dict(status='passed',old_files_unchanged=same,
        changed_entries=changed,prior_bound_files=len(prior['hashes']),old_archives=archives),indent=2))
    sources={
        'sources/kim-scarinci-2112.13329.pdf':'https://arxiv.org/pdf/2112.13329v2',
        'sources/ramond-distributions-2015.pdf':'https://www.math.u-psud.fr/~ramond/docs/cours/distributionsedp.pdf'}
    (base/'source-record.json').write_text(json.dumps(dict(sources=[dict(file=p,url=u,sha256=sha(root/p)) for p,u in sources.items()],
        prior_evidence={priorpath.relative_to(root).as_posix():sha(priorpath)},
        reading_scope='KS v2 scalar product and pole divisor; Ramond version 2.06, Example 4.4.5 and Proposition 4.4.6, printed p87 / PDF p88 rendered. Standard distribution identity, not a new general theorem.'),indent=2))
    replay=audit(root);assert replay['status']=='passed' and replay['checks']==62,replay
    paths=['code/'+p for p in ['weak_continuation.py','run_weak_continuation.py','positive_weak_form.py',
           'check_weak_continuation.py','audit_weak_continuation.py','build_weak_continuation.ps1','record_weak_continuation.py']]
    paths+=['manuscript/weak-continuation/main.tex','output/pdf/de-sitter-weak-continuation.pdf',
            'notes/stage4-weak-continuation.md','literature/stage4-weak-continuation-sources.md']
    paths+=[p.relative_to(root).as_posix() for p in sorted(base.iterdir()) if p.is_file()]
    evidence=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),status='diagnostic stage bound; project audit follows',
        scientific_status='literal local weak pairing is not holomorphic; bounded pentagon not tested',pages=5,
        runs=13,positive_state_records=3,diagnostic_checks=34,saved_arithmetic_replay_checks=62,
        compilers=['pdfLaTeX twice exit 0','desktop success'],
        visual_review=dict(pages=list(range(1,6)),longest_edge_pixels=1100,status='all inspected'),
        hashes={p:sha(root/p) for p in paths})
    (base/'evidence.json').write_text(json.dumps(evidence,indent=2))
    print(json.dumps(dict(bound_files=len(paths),old_files_unchanged=same,old_archives=len(archives),scientific_status=evidence['scientific_status'])))
if __name__=='__main__':main()
