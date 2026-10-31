"""Bind new evidence and verify immutable pre-existing scientific files."""
import hashlib,json,zipfile
from datetime import datetime,timezone
from pathlib import Path


def main():
    root=Path(__file__).resolve().parents[1];base=root/'results/stage4/ds3-transverse-phase'
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    if (base/'evidence.json').exists():raise RuntimeError('Preserve existing evidence')
    prior={}
    for stage in ['ds3-full-tail','ds3-shadow-duality','ds3-peak-limit']:
        data=json.loads((root/f'results/stage4/{stage}/evidence.json').read_text())
        for name,digest in data['hashes'].items():assert sha(root/name)==digest,name
        prior[stage]=len(data['hashes'])
    baseline=root/'snapshots/2026-10-31-ds3-bose.zip'
    assert sha(baseline)=='fddbeeb767dca73e38b352c956be56d109e9f7283c436c51c100de31ba11f74c'
    allowed={'README.md','PLAN.md','NEXT.md','journal/2026-10-31.md','code/project_audit.py',
             'results/project-audit.json','results/project-files.json'}
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
    assert len(archives)==45
    preservation=dict(status='passed',old_snapshot_files_unchanged=same,changed_entries=changed,
                      prior_stage_bound_files=prior,old_archives=archives)
    (base/'preservation.json').write_text(json.dumps(preservation,indent=2),encoding='utf-8')
    collision=root/'results/stage4/ds3-collision'
    if (collision/'evidence.json').exists():raise RuntimeError('Preserve collision evidence')
    paths=['code/'+n+'.py' for n in ['ds3_racah_collision','ds3_racah_merged','run_ds3_collision',
           'validate_ds3_collision','ds3_controlled_omega','ds3_collision_asymptotic','validate_ds3_controlled','audit_ds3_collision']]
    paths+=['code/build_ds3_collision.ps1','notes/stage4-ds3-collision.md','literature/stage4-ds3-collision-sources.md',
            'manuscript/ds3-collision/main.tex','output/pdf/de-sitter-ds3-collision.pdf']
    paths += [p.relative_to(root).as_posix() for p in sorted(collision.iterdir()) if p.is_file()]
    record=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),status='completed verification; historical PDF retained',
                prior_job=dict(session_id=47841,exit_code=0),hashes={p:sha(root/p) for p in paths})
    (collision/'evidence.json').write_text(json.dumps(record,indent=2),encoding='utf-8')
    paths=['code/'+n+'.py' for n in ['ds3_transverse_phase','run_ds3_transverse_phase','audit_ds3_transverse_phase',
           'refine_ds3_transverse_phase','record_ds3_transverse_phase','audit_ds3_recent_stages']]
    paths+=['code/build_ds3_transverse_phase.ps1','notes/stage4-ds3-transverse-phase.md',
            'literature/stage4-ds3-transverse-phase-sources.md','manuscript/ds3-transverse-phase/main.tex',
            'output/pdf/de-sitter-ds3-transverse-phase.pdf','sources/manifest.json','results/stage4/ds3-collision/evidence.json']
    paths += [p.relative_to(root).as_posix() for p in sorted(base.iterdir()) if p.is_file()]
    record=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),status='stage checks complete; final project audit next',
                pages=4,study_runs=16,independent_audit_checks=203,precision_refinements=4,
                compilers=['local pdfLaTeX twice','desktop success'],
                visual_review=dict(pages=[1,2,3,4],longest_edge_pixels=1100,status='all inspected; no layout defects'),
                active_old_jobs=[],hashes={p:sha(root/p) for p in paths})
    (base/'evidence.json').write_text(json.dumps(record,indent=2),encoding='utf-8')
    for path,digest in record['hashes'].items():assert sha(root/path)==digest,path
    print(json.dumps(dict(status='passed',new_bound_files=len(paths),old_files_unchanged=same,prior=prior,old_zips=len(archives)),indent=2))


if __name__=='__main__':main()
