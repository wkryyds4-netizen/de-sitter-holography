"""Bind a limited domain/literature diagnostic, preserving earlier evidence."""
from pathlib import Path
import json,hashlib,zipfile
from datetime import datetime,timezone
from audit_flip_domain import audit

def main():
    root=Path(__file__).resolve().parents[1];base=root/'results/stage4/flip-domain'
    if (base/'evidence.json').exists():raise RuntimeError('Already bound')
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    priorpath=root/'results/stage4/unitary-flip/evidence.json'
    prior=json.loads(priorpath.read_text())
    for p,h in prior['hashes'].items():assert sha(root/p)==h,p
    known={r['name']:r['sha256'] for r in json.loads((root/'results/stage4/unitary-flip/preservation.json').read_text())['old_archives']}
    archives=[]
    for p in sorted((root/'snapshots').glob('*.zip')):
        h=sha(p);side=p.with_suffix('.sha256');assert side.exists() or p.name in known,p
        if side.exists():assert side.read_text().split()[0]==h,p
        if p.name in known:assert h==known[p.name],p
        archives.append(dict(name=p.name,sha256=h))
    assert len(archives)==54 and len(known)==53
    same=0;changed=[]
    allowed={'README.md','PLAN.md','NEXT.md','journal/2026-11-01.md','code/project_audit.py','results/project-audit.json','results/project-files.json'}
    with zipfile.ZipFile(root/'snapshots/2026-11-01-unitary-flip.zip') as z:
        for p in z.namelist():
            if p.endswith('/'):continue
            if (root/p).is_file() and sha(root/p)==hashlib.sha256(z.read(p)).hexdigest():same+=1
            else:assert p in allowed,p;changed.append(p)
    (base/'preservation.json').write_text(json.dumps(dict(status='passed',old_files_unchanged=same,
        changed_entries=changed,prior_bound_files=len(prior['hashes']),old_archives=archives),indent=2))
    sources={
        'sources/kim-scarinci-2112.13329.pdf':'https://arxiv.org/pdf/2112.13329v2',
        'sources/kim-trilogy-2308.02934v2.html':'https://arxiv.org/html/2308.02934v2',
        'sources/fock-goncharov-math-0702397v6.pdf':'https://arxiv.org/pdf/math/0702397v6',
        'sources/goncharov-0706.4054v2.pdf':'https://arxiv.org/pdf/0706.4054v2',
        'sources/kim-pseudo-kahler-2112.13962v2.pdf':'https://arxiv.org/pdf/2112.13962v2',
        'sources/pusz-soltan-math-0604613v1.pdf':'https://arxiv.org/pdf/math/0604613v1'}
    (base/'source-record.json').write_text(json.dumps(dict(sources=[dict(file=p,url=u,sha256=sha(root/p)) for p,u in sources.items()],
        prior_evidence={priorpath.relative_to(root).as_posix():sha(priorpath)},
        reading_scope='Pinned local partial reading. FG p38, Kim pseudo-Kahler p21, Pusz-Soltan p3 rendered; no full independent theorem verification.'),indent=2))
    replay=audit(root);assert replay['status']=='passed' and replay['checks']==322,replay
    paths=['code/'+p for p in ['flip_domain.py','run_flip_domain.py','check_flip_domain.py','refine_flip_domain.py',
           'check_flip_regularity.py','audit_flip_domain.py','build_flip_domain.ps1','record_flip_domain.py']]
    paths+=['manuscript/flip-domain/main.tex','output/pdf/de-sitter-flip-domain.pdf','notes/stage4-flip-domain.md',
            'literature/stage4-flip-domain-sources.md','literature/flip-proof-dependency-map.md']
    paths+=[p.relative_to(root).as_posix() for p in sorted(base.iterdir()) if p.is_file()]
    evidence=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),status='diagnostic stage bound; project audit follows',
        scientific_status='bounded pentagon unresolved; quadrature not certified',pages=6,
        runs=11,initial_checks=dict(total=22,passed=20,failed=2),refinement_checks=5,regularity_checks=31,
        saved_vector_replay_checks=322,compilers=['pdfLaTeX twice exit 0','desktop success'],
        visual_review=dict(pages=list(range(1,7)),longest_edge_pixels=1100,status='all inspected'),
        hashes={p:sha(root/p) for p in paths})
    (base/'evidence.json').write_text(json.dumps(evidence,indent=2))
    print(json.dumps(dict(bound_files=len(paths),old_files_unchanged=same,old_archives=len(archives),scientific_status=evidence['scientific_status'])))
if __name__=='__main__':main()
