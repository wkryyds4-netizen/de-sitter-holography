"""Bind this stage and verify preservation; never rewrite an earlier record."""
import hashlib,json,zipfile
from datetime import datetime,timezone
from pathlib import Path
from audit_principal_eikonal import audit

def main():
    root=Path(__file__).resolve().parents[1];base=root/'results/stage4/principal-eikonal'
    if (base/'evidence.json').exists():raise RuntimeError('Already frozen')
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    prior=json.loads((root/'results/stage4/gaussian-eikonal/evidence.json').read_text())
    for p,h in prior['hashes'].items():assert sha(root/p)==h,p
    archives=[]
    previous_archives=json.loads((root/'results/stage4/gaussian-eikonal/preservation.json').read_text())['old_archives']
    known={r['name']:r['sha256'] for r in previous_archives}
    for path in sorted((root/'snapshots').glob('*.zip')):
        digest=sha(path);side=path.with_suffix('.sha256')
        assert side.exists() or path.name in known,path
        if side.exists():assert digest==side.read_text().split()[0],path
        if path.name in known:assert digest==known[path.name],path
        archives.append(dict(name=path.name,sha256=digest))
    assert len(archives)==51 and len(known)==50
    same=0;changed=[]
    allowed={'README.md','PLAN.md','NEXT.md','journal/2026-10-31.md','code/project_audit.py','results/project-audit.json','results/project-files.json'}
    with zipfile.ZipFile(root/'snapshots/2026-10-31-gaussian-eikonal.zip') as z:
        for name in z.namelist():
            if name.endswith('/'):continue
            if (root/name).is_file() and sha(root/name)==hashlib.sha256(z.read(name)).hexdigest():same+=1
            else:assert name in allowed,name;changed.append(name)
    (base/'preservation.json').write_text(json.dumps(dict(status='passed',old_files_unchanged=same,
        changed_entries=changed,prior_bound_files=len(prior['hashes']),old_archives=archives),indent=2))
    sourcefiles={
        'sources/chen-stanford-tang-yang-2607.14042v1.html':'https://arxiv.org/html/2607.14042v1',
        'sources/milekhin-zero-2607.13137v1.html':'https://arxiv.org/html/2607.13137v1',
        'sources/bonifacio-edge-2610.28719v1.pdf':'https://arxiv.org/pdf/2610.28719v1'}
    sources=[dict(file=p,url=u,sha256=sha(root/p)) for p,u in sourcefiles.items()]
    (base/'source-record.json').write_text(json.dumps(dict(status='local caches checked; full texts excluded from Git',sources=sources,
        online_formulas=['https://dlmf.nist.gov/10.32.E9','https://dlmf.nist.gov/10.43.E19']),indent=2))
    a=audit(root);assert a['status']=='passed' and a['checks']==727,a
    code=['principal_eikonal.py','run_principal_eikonal.py','check_principal_eikonal.py','refine_principal_eikonal.py',
          'late_principal_eikonal.py','audit_principal_eikonal.py','record_principal_eikonal.py','build_principal_eikonal.ps1']
    paths=['code/'+p for p in code]+['manuscript/principal-eikonal/main.tex','output/pdf/de-sitter-principal-eikonal.pdf',
          'notes/stage4-principal-eikonal.md','literature/stage4-principal-eikonal-sources.md']
    paths += [p.relative_to(root).as_posix() for p in sorted(base.iterdir()) if p.is_file()]
    result=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),status='stage verified; project audit follows',pages=6,
        full_scrambled_runs=144,total_sampled_points=184549376,diagnostic_counts=[60,14,727],
        compilers=['pdfLaTeX twice exit 0','desktop success'],
        visual_review=dict(pages=list(range(1,7)),longest_edge_pixels=1100,status='all pages inspected'),
        scope='Specified real-principal eikonal lateral kernel, Gaussian moment spectrum and signed full model integrals. Empirical quadrature errors; late u8 sign unresolved. Not full finite-G dynamics, physical state identification or microscopic dual.',
        hashes={p:sha(root/p) for p in paths})
    (base/'evidence.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(dict(status='passed',bound_files=len(paths),old_files_unchanged=same,old_archives=len(archives))))

if __name__=='__main__':main()
