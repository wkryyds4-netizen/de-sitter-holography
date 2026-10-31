"""Freeze new evidence without modifying prior scientific results."""
import hashlib,json,zipfile
from datetime import datetime,timezone
from pathlib import Path


def main():
    root=Path(__file__).resolve().parents[1];base=root/'results/stage4/ds3-physical-local'
    if (base/'evidence.json').exists():raise RuntimeError('Evidence already frozen')
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    prior=json.loads((root/'results/stage4/ds3-uniform-local/evidence.json').read_text())
    for name,digest in prior['hashes'].items():assert sha(root/name)==digest,name
    snap=root/'snapshots/2026-10-31-ds3-uniform-local.zip'
    assert sha(snap)==snap.with_suffix('.sha256').read_text().split()[0]
    allowed={'README.md','PLAN.md','NEXT.md','REPOSITORY.md','journal/2026-10-31.md',
             'code/project_audit.py','results/project-audit.json','results/project-files.json'}
    same=0;changed=[]
    with zipfile.ZipFile(snap) as z:
        for name in z.namelist():
            if name.endswith('/'):continue
            if (root/name).is_file() and sha(root/name)==hashlib.sha256(z.read(name)).hexdigest():same+=1
            else:assert name in allowed,name;changed.append(name)
    archives=[]
    for path in sorted((root/'snapshots').glob('*.zip')):
        digest=sha(path);side=path.with_suffix('.sha256')
        if side.exists():assert side.read_text().split()[0]==digest,path
        archives.append(dict(name=path.name,sha256=digest))
    (base/'preservation.json').write_text(json.dumps(dict(status='passed',old_files_unchanged=same,
        changed_entries=changed,prior_bound_files=len(prior['hashes']),old_archives=archives),indent=2),encoding='utf-8')
    paths=['code/'+s for s in ['ds3_physical_local.py','check_ds3_physical_local.py',
           'audit_ds3_physical_local.py','record_ds3_physical_local.py','build_ds3_physical_local.ps1']]
    paths+=['manuscript/ds3-physical-local/main.tex','output/pdf/de-sitter-ds3-physical-local.pdf',
            'notes/stage4-ds3-physical-local.md','literature/stage4-ds3-physical-local-sources.md']
    paths+=[p.relative_to(root).as_posix() for p in sorted(base.iterdir()) if p.is_file()]
    result=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),status='stage verified; project audit follows',
                pages=4,diagnostic_checks=97,independent_decimal_checks=165,reused_full_kernel_nodes=144,
                new_full_kernel_evaluations=0,compilers=['local pdfLaTeX twice, exit 0','desktop success'],
                visual_review=dict(pages=[1,2,3,4],longest_edge_pixels=1100,status='all pages inspected'),
                scope='Compact physical weights and cutoff pairing; full Bose/free limit remains open',
                hashes={p:sha(root/p) for p in paths})
    (base/'evidence.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(dict(status='passed',bound_files=len(paths),old_files_unchanged=same,old_archives=len(archives))))


if __name__=='__main__':main()
