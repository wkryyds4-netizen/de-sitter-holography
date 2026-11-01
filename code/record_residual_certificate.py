"""Bind completed residual work and preserve the previous 58 checkpoints."""
from pathlib import Path
import json,hashlib,zipfile
from datetime import datetime,timezone
from audit_tensor_chain import audit
from audit_chirp_residual import audit as chirp_audit

def main():
    root=Path(__file__).resolve().parents[1];base=root/'results/stage4/residual-certificate'
    if (base/'evidence.json').exists():raise RuntimeError('Already bound')
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    priorpath=root/'results/stage4/heat-matrix/evidence.json';prior=json.loads(priorpath.read_text())
    for p,h in prior['hashes'].items():assert sha(root/p)==h,p
    known={r['name']:r['sha256'] for r in json.loads((root/'results/stage4/heat-matrix/preservation.json').read_text())['old_archives']}
    archives=[]
    for p in sorted((root/'snapshots').glob('*.zip')):
        digest=sha(p);side=p.with_suffix('.sha256')
        assert side.exists() or p.name in known,p
        if side.exists():assert side.read_text().split()[0]==digest,p
        if p.name in known:assert known[p.name]==digest,p
        archives.append(dict(name=p.name,sha256=digest))
    assert len(archives)==58 and len(known)==57
    same=0;changed=[]
    allowed={'README.md','PLAN.md','NEXT.md','journal/2026-11-01.md','code/project_audit.py','results/project-audit.json','results/project-files.json'}
    with zipfile.ZipFile(root/'snapshots/2026-11-01-heat-matrix.zip') as z:
        for p in z.namelist():
            if p.endswith('/'):continue
            if (root/p).is_file() and sha(root/p)==hashlib.sha256(z.read(p)).hexdigest():same+=1
            else:assert p in allowed,p;changed.append(p)
    with (base/'preservation.json').open('x') as f:json.dump(dict(status='passed',old_files_unchanged=same,changed_entries=changed,prior_bound_files=len(prior['hashes']),old_archives=archives),f,indent=2)
    sources={'sources/kim-scarinci-2112.13329.pdf':'https://arxiv.org/pdf/2112.13329v2'}
    online=['https://arxiv.org/html/2412.07562v2','https://arxiv.org/pdf/1409.1208v1','https://arxiv.org/html/1802.07942v1','https://dlmf.nist.gov/18.18.E28','https://dlmf.nist.gov/18.9','https://arxiv.org/html/2308.02934v2','https://python-flint.readthedocs.io/en/latest/acb.html']
    with (base/'source-record.json').open('x') as f:json.dump(dict(sources=[dict(file=p,url=u,sha256=sha(root/p)) for p,u in sources.items()],online_references=online,prior_evidence={priorpath.relative_to(root).as_posix():sha(priorpath)},reading_scope='KS sections 6.3 and Proposition 6.15; online localized formula/assumption checks, not full-paper verification or exhaustive novelty search.'),f,indent=2)
    chain=audit(root,replay_bounds=False);chirp=chirp_audit()
    assert chain['passed'] and chirp['status']=='passed'
    names=['chirp_residual_certificate.py','chirp_residual_certificate_initial.py','audit_chirp_residual.py','residual_overlap_rectangle.py','residual_full_overlap.py','residual_gaussian_overlap.py','residual_gaussian_overlap_bounded.py','gaussian_outer_bound.py','tensor_overlap.py','run_tensor_gaussian.py','check_tensor_overlap.py','audit_tensor_controls.py','audit_tensor_gaussian.py','audit_tensor_gaussian_refined.py','tensor_pair_overlap.py','check_tensor_pair.py','run_tensor_pair.py','audit_tensor_chain.py','build_residual_certificate.ps1','record_residual_certificate.py','requirements-residual.txt']
    names.append('check_actual_pair_rectangle.py')
    paths=['code/'+p for p in names]+['manuscript/residual-certificate/main.tex','output/pdf/de-sitter-residual-certificate.pdf','notes/stage4-residual-certificate.md','literature/stage4-residual-certificate-sources.md']
    paths+=[p.relative_to(root).as_posix() for p in sorted(base.iterdir()) if p.is_file()]
    review=json.loads((base/'final-build.json').read_text())
    record=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),status='stage bound; global audit follows',scientific_status=chain['scientific_status'],scope=chain['scope'],pages=review['pages'],compilers=['pdfLaTeX twice exit 0'],visual_review=review['visual_review'],word_norm_lower_bound=chain['word_norm_lower_bound'],hashes={p:sha(root/p) for p in paths})
    with (base/'evidence.json').open('x') as f:json.dump(record,f,indent=2)
    print(json.dumps(dict(bound_files=len(paths),old_files_unchanged=same,archives=len(archives),scientific_status=record['scientific_status'])))

if __name__=='__main__':main()
