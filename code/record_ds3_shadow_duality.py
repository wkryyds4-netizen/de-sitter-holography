"""Bind the new stage and verify all previously frozen scientific bytes."""
import hashlib,json,zipfile
from datetime import datetime,timezone
from pathlib import Path


def main():
    root=Path(__file__).resolve().parents[1];base=root/'results/stage4/ds3-shadow-duality';dest=base/'evidence.json'
    if dest.exists():raise RuntimeError('Preserve existing evidence')
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    old=json.loads((root/'results/stage4/ds3-full-tail/evidence.json').read_text())
    for name,digest in old['hashes'].items():assert sha(root/name)==digest,name
    baseline=root/'snapshots/2026-10-31-ds3-bose.zip'
    assert sha(baseline)=='fddbeeb767dca73e38b352c956be56d109e9f7283c436c51c100de31ba11f74c'
    allowed={'README.md','PLAN.md','NEXT.md','journal/2026-10-31.md','results/project-audit.json','results/project-files.json'}
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
                      previous_unfrozen_stage_hashes_checked=len(old['hashes']),old_archives=archives)
    (base/'preservation.json').write_text(json.dumps(preservation,indent=2),encoding='utf-8')
    paths=['code/'+n+'.py' for n in ['ds3_shadow_duality','validate_ds3_shadow_duality','run_ds3_shadow_duality','refine_ds3_shadow_identity','audit_ds3_shadow_duality','record_ds3_shadow_duality']]
    paths+=['code/build_ds3_shadow_duality.ps1','notes/stage4-ds3-shadow-duality.md','literature/stage4-ds3-shadow-duality-sources.md','manuscript/ds3-shadow-duality/main.tex','output/pdf/de-sitter-ds3-shadow-duality.pdf','sources/manifest.json']
    paths += [p.relative_to(root).as_posix() for p in sorted(base.iterdir()) if p.is_file()]
    evidence=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),status='stage checks completed; project integration and freeze pending',
                  pages=4,initial_checks=dict(total=33,passed=31,failed=2,preserved=True),fresh_refinements=4,independent_audit_checks=365,
                  study_runs=18,compilers=['local pdfLaTeX twice','desktop success'],visual_review=dict(pages=[1,2,3,4],longest_edge_pixels=1100,status='all inspected, no layout defects'),
                  pending_prior_job=dict(session_id=47841,counted_as_passed=False),hashes={p:sha(root/p) for p in paths})
    dest.write_text(json.dumps(evidence,indent=2),encoding='utf-8')
    for p,digest in evidence['hashes'].items():assert sha(root/p)==digest,p
    print(json.dumps(dict(status='passed',evidence_files=len(paths),old_files_unchanged=same,prior_stage_bound_files=len(old['hashes']),zip_count=len(archives)),indent=2))


if __name__=='__main__':main()
