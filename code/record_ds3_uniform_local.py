"""Record the local-limit stage, preserving prior scientific bytes and archives."""
import hashlib,json,zipfile
from datetime import datetime,timezone
from pathlib import Path


def main():
    root=Path(__file__).resolve().parents[1];base=root/'results/stage4/ds3-uniform-local'
    if (base/'evidence.json').exists():raise RuntimeError('Preserve existing evidence')
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    prior={}
    for stage in ['ds3-collision','ds3-full-tail','ds3-shadow-duality','ds3-peak-limit','ds3-transverse-phase']:
        data=json.loads((root/f'results/stage4/{stage}/evidence.json').read_text())
        for name,digest in data['hashes'].items():assert sha(root/name)==digest,name
        prior[stage]=len(data['hashes'])
    baseline=root/'snapshots/2026-10-31-ds3-transverse-phase.zip'
    assert sha(baseline)=='da0b5a4e7126f5bc6de432d66755cd0048043911e388b846126af08d85040d11'
    allowed={'README.md','PLAN.md','NEXT.md','journal/2026-10-31.md','code/project_audit.py',
             'results/project-audit.json','results/project-files.json','.gitignore','AGENTS.md'}
    same=0;changed=[]
    with zipfile.ZipFile(baseline) as archive:
        for name in archive.namelist():
            if name.endswith('/'):continue
            if (root/name).is_file() and sha(root/name)==hashlib.sha256(archive.read(name)).hexdigest():same+=1
            else:assert name in allowed,name;changed.append(name)
    archives=[]
    for p in sorted((root/'snapshots').glob('*.zip')):
        digest=sha(p);side=p.with_suffix('.sha256')
        if side.exists():assert side.read_text().split()[0]==digest,p
        archives.append(dict(file=p.name,sha256=digest,sidecar_present=side.exists()))
    assert len(archives)==46
    preservation=dict(status='passed',old_snapshot_files_unchanged=same,changed_entries=changed,
                      prior_bound_files=prior,old_archives=archives,
                      context_changes='AGENTS refreshed by user; .gitignore updated for explicitly requested Git publication.')
    (base/'preservation.json').write_text(json.dumps(preservation,indent=2),encoding='utf-8')
    paths=['code/'+n+'.py' for n in ['ds3_uniform_local','run_ds3_uniform_local','validate_ds3_uniform_local',
           'integrate_ds3_uniform_local','audit_ds3_uniform_local','record_ds3_uniform_local']]
    paths+=['code/build_ds3_uniform_local.ps1','notes/stage4-ds3-uniform-local.md',
            'literature/stage4-ds3-uniform-local-sources.md','manuscript/ds3-uniform-local/main.tex',
            'output/pdf/de-sitter-ds3-uniform-local.pdf','sources/manifest.json']
    paths+=[p.relative_to(root).as_posix() for p in sorted(base.iterdir()) if p.is_file()]
    result=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),status='stage verification complete; project audit follows',
                pages=5,study_points=21,independent_representation_checks=11,quadrature_nodes=144,audit_checks=446,
                compilers=['local pdfLaTeX twice','desktop success'],
                visual_review=dict(pages=[1,2,3,4,5],longest_edge_pixels=1100,status='all inspected; final page 2 condition update rechecked'),
                hashes={p:sha(root/p) for p in paths})
    (base/'evidence.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    for name,digest in result['hashes'].items():assert sha(root/name)==digest,name
    print(json.dumps(dict(status='passed',bound_files=len(paths),old_files_unchanged=same,old_zips=len(archives),prior=prior),indent=2))


if __name__=='__main__':main()
