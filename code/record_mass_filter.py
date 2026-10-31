"""Bind this stage and verify preservation; never rewrite an earlier record."""
import hashlib,json,zipfile
from datetime import datetime,timezone
from pathlib import Path
from audit_mass_filter import audit

def main():
    root=Path(__file__).resolve().parents[1];base=root/'results/stage4/mass-filter'
    if (base/'evidence.json').exists():raise RuntimeError('Already frozen')
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    prior=json.loads((root/'results/stage4/principal-eikonal/evidence.json').read_text())
    for p,h in prior['hashes'].items():assert sha(root/p)==h,p
    archives=[]
    previous_archives=json.loads((root/'results/stage4/principal-eikonal/preservation.json').read_text())['old_archives']
    known={r['name']:r['sha256'] for r in previous_archives}
    for path in sorted((root/'snapshots').glob('*.zip')):
        digest=sha(path);side=path.with_suffix('.sha256')
        assert side.exists() or path.name in known,path
        if side.exists():assert digest==side.read_text().split()[0],path
        if path.name in known:assert digest==known[path.name],path
        archives.append(dict(name=path.name,sha256=digest))
    assert len(archives)==52 and len(known)==51
    same=0;changed=[]
    allowed={'README.md','PLAN.md','NEXT.md','journal/2026-10-31.md','code/project_audit.py','results/project-audit.json','results/project-files.json'}
    with zipfile.ZipFile(root/'snapshots/2026-10-31-principal-eikonal.zip') as z:
        for name in z.namelist():
            if name.endswith('/'):continue
            if (root/name).is_file() and sha(root/name)==hashlib.sha256(z.read(name)).hexdigest():same+=1
            else:assert name in allowed,name;changed.append(name)
    (base/'preservation.json').write_text(json.dumps(dict(status='passed',old_files_unchanged=same,
        changed_entries=changed,prior_bound_files=len(prior['hashes']),old_archives=archives),indent=2))
    sourcefiles={
        'sources/chen-stanford-tang-yang-2607.14042v1.html':'https://arxiv.org/html/2607.14042v1',
        'sources/chandra-liska-post-2610.30247v1.html':'https://arxiv.org/html/2610.30247v1'}
    sources=[dict(file=p,url=u,sha256=sha(root/p)) for p,u in sourcefiles.items()]
    priorfiles=['results/stage4/principal-eikonal/evidence.json','results/stage4/recoil-spectral/evidence.json']
    (base/'source-record.json').write_text(json.dumps(dict(status='local caches checked; full texts excluded from Git',sources=sources,
        prior_evidence={p:sha(root/p) for p in priorfiles}),indent=2))
    a=audit(root);assert a['status']=='passed' and a['checks']==444,a
    code=['mass_filter.py','run_mass_filter.py','check_mass_filter.py','packet_mass_filter.py',
          'audit_mass_filter.py','record_mass_filter.py','build_mass_filter.ps1']
    paths=['code/'+p for p in code]+['manuscript/mass-filter/main.tex','output/pdf/de-sitter-mass-filter.pdf',
          'notes/stage4-mass-filter.md','literature/stage4-mass-filter-sources.md']
    paths += [p.relative_to(root).as_posix() for p in sorted(base.iterdir()) if p.is_file()]
    result=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),status='stage verified; project audit follows',pages=5,
        full_scrambled_runs=80,total_sampled_points=104857600,diagnostic_counts=[23,444],
        compilers=['pdfLaTeX twice exit 0','desktop success'],
        visual_review=dict(pages=list(range(1,6)),longest_edge_pixels=1100,status='all pages inspected'),
        scope='Exact Gaussian filter transform and specified eikonal comparison models. Finite-kappa mass map is an additional assumption; not an established physical trace, full finite-G dynamics or microscopic dual. Errors are empirical.',
        hashes={p:sha(root/p) for p in paths})
    (base/'evidence.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(dict(status='passed',bound_files=len(paths),old_files_unchanged=same,old_archives=len(archives))))

if __name__=='__main__':main()
