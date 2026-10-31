"""Bind Gaussian eikonal evidence without rewriting earlier scientific records."""
import hashlib,json,zipfile
from datetime import datetime,timezone
from pathlib import Path

def main():
    root=Path(__file__).resolve().parents[1];base=root/'results/stage4/gaussian-eikonal'
    if (base/'evidence.json').exists():raise RuntimeError('Evidence already frozen')
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    prior=json.loads((root/'results/stage4/ds3-fock-vacuum/evidence.json').read_text())
    for name,digest in prior['hashes'].items():assert sha(root/name)==digest,name
    snap=root/'snapshots/2026-10-31-ds3-fock-vacuum.zip'
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
    assert len(archives)==50,'Expected 50 prior snapshots'
    (base/'preservation.json').write_text(json.dumps(dict(status='passed',old_files_unchanged=same,
        changed_entries=changed,prior_bound_files=len(prior['hashes']),old_archives=archives),indent=2),encoding='utf-8')
    paths=['code/'+name for name in ['gaussian_eikonal.py','run_gaussian_eikonal.py',
        'check_gaussian_eikonal.py','refine_gaussian_eikonal.py','record_gaussian_eikonal.py','build_gaussian_eikonal.ps1']]
    paths+=['manuscript/gaussian-eikonal/main.tex','output/pdf/de-sitter-gaussian-eikonal.pdf',
        'notes/stage4-gaussian-eikonal.md','literature/stage4-gaussian-eikonal-sources.md']
    paths+=[p.relative_to(root).as_posix() for p in sorted(base.iterdir()) if p.is_file()]
    for filename,count in [('diagnostics.json',52),('refinement.json',11)]:
        d=json.loads((base/filename).read_text());assert d['status']=='passed' and len(d['checks'])==count and all(x['passed'] for x in d['checks'])
    result=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),status='stage verified; project audit follows',
        pages=6,eikonal_parameter_points=8,diagnostic_counts=[52,11],
        original_four_coordinate_controls=2,original_momentum_integrals=3,
        compilers=['local pdfLaTeX twice, exit 0','desktop success'],
        visual_review=dict(pages=list(range(1,7)),longest_edge_pixels=1100,status='all pages inspected; final bibliography correction rechecked on page 6'),
        scope='Four independent Gaussian smearings of the specified two-conformal-probe eikonal contours, normalization and surviving imaginary sign. Not full finite-G correlators, principal-series Fock continuation or a microscopic dual. Gaussian tails do not bound omitted gravity terms.',
        hashes={p:sha(root/p) for p in paths})
    (base/'evidence.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(dict(status='passed',bound_files=len(paths),old_files_unchanged=same,old_archives=len(archives))))

if __name__=='__main__':main()
