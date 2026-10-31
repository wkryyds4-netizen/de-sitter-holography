"""Bind regulated Fock evidence without rewriting earlier scientific records."""
import hashlib,json,zipfile
from datetime import datetime,timezone
from pathlib import Path

def main():
    root=Path(__file__).resolve().parents[1];base=root/'results/stage4/ds3-fock-vacuum'
    if (base/'evidence.json').exists():raise RuntimeError('Evidence already frozen')
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    prior=json.loads((root/'results/stage4/ds3-compact-global/evidence.json').read_text())
    for name,digest in prior['hashes'].items():assert sha(root/name)==digest,name
    snap=root/'snapshots/2026-10-31-ds3-compact-global.zip'
    assert sha(snap)==snap.with_suffix('.sha256').read_text().split()[0]
    allowed={'README.md','PLAN.md','NEXT.md','journal/2026-10-31.md',
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
    assert len(archives)==49,'Expected 49 prior snapshots'
    (base/'preservation.json').write_text(json.dumps(dict(status='passed',old_files_unchanged=same,
        changed_entries=changed,prior_bound_files=len(prior['hashes']),old_archives=archives),indent=2),encoding='utf-8')
    paths=['code/'+s for s in ['ds3_fock_vacuum.py','check_ds3_fock_vacuum.py',
        'recover_ds3_fock_vacuum.py','audit_ds3_fock_vacuum.py','audit_ds3_fock_vacuum_refined.py',
        'ds3_gaussian_field.py','check_ds3_gaussian_field.py','audit_ds3_gaussian_field.py',
        'record_ds3_fock_vacuum.py','build_ds3_fock_vacuum.ps1']]
    paths+=['manuscript/ds3-fock-vacuum/main.tex','output/pdf/de-sitter-ds3-fock-vacuum.pdf',
        'notes/stage4-ds3-fock-vacuum.md','literature/stage4-ds3-fock-vacuum-sources.md']
    paths+=[p.relative_to(root).as_posix() for p in sorted(base.iterdir()) if p.is_file()]
    for filename,count in [('audit-refined.json',154),('gaussian-audit.json',115)]:
        d=json.loads((base/filename).read_text());assert d['status']=='passed' and d['checks']==count and not d['errors']
    result=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),status='stage verified; project audit follows',
        pages=9,balanced_return_pairs=20,gaussian_return_pairs=16,diagnostic_counts=[36,24],
        independent_check_counts=[154,115],full_crossed_kernel_evaluations=0,full_fourfield_quadratures=0,
        compilers=['local pdfLaTeX twice, exit 0','desktop success'],
        visual_review=dict(pages=list(range(1,10)),longest_edge_pixels=1100,status='all pages inspected'),
        failures_retained=['NumPy Boolean JSON serialization after raw quadrature completion',
                           'SymPy structural equality; exact difference repair retained separately'],
        scope='Self-checked regulated G=0 full multiplicity Fock completion, weak descent and vacuum fourpoint limit. Gaussian thermal spectrum differs from balanced completion. No unregulated field, finite-gravity dynamics, factor type or full holographic dual claim.',
        hashes={p:sha(root/p) for p in paths})
    (base/'evidence.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(dict(status='passed',bound_files=len(paths),old_files_unchanged=same,old_archives=len(archives))))

if __name__=='__main__':main()
