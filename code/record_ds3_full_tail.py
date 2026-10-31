"""Save this unfrozen stage's evidence and verify preservation of the last snapshot."""
import hashlib
import json
import zipfile
from datetime import datetime, timezone
from pathlib import Path


def main():
    root=Path(__file__).resolve().parents[1];base=root/'results/stage4/ds3-full-tail'
    target=base/'evidence.json'
    if target.exists():raise RuntimeError('Preserve prior evidence')
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    baseline=root/'snapshots/2026-10-31-ds3-bose.zip'
    assert sha(baseline)=='fddbeeb767dca73e38b352c956be56d109e9f7283c436c51c100de31ba11f74c'
    allowed={'README.md','PLAN.md','NEXT.md','journal/2026-10-31.md','results/project-audit.json','results/project-files.json'}
    unchanged=[];changed=[];unexpected=[]
    with zipfile.ZipFile(baseline) as archive:
        for name in archive.namelist():
            if name.endswith('/'):continue
            p=root/name
            same=p.is_file() and sha(p)==hashlib.sha256(archive.read(name)).hexdigest()
            if same:unchanged.append(name)
            elif name in allowed:changed.append(name)
            else:unexpected.append(name)
    archives=[]
    for path in sorted((root/'snapshots').glob('*.zip')):
        digest=sha(path);side=path.with_suffix('.sha256')
        matched=digest==side.read_text().split()[0] if side.exists() else None
        if matched is False:unexpected.append(path.relative_to(root).as_posix())
        archives.append(dict(file=path.relative_to(root).as_posix(),sha256=digest,sidecar_matches=matched))
    assert len(archives)==45 and not unexpected,(len(archives),unexpected)
    preservation=dict(status='passed',baseline=baseline.relative_to(root).as_posix(),unchanged_files=len(unchanged),
                      changed_entry_files=changed,unexpected_changes=unexpected,archives=archives)
    (base/'preservation.json').write_text(json.dumps(preservation,indent=2),encoding='utf-8')
    paths=[f'code/{n}.py' for n in ['ds3_full_tail_bound','run_ds3_full_tail','analyze_ds3_shadow_amplification','audit_ds3_full_tail','record_ds3_full_tail']]
    paths+=['code/build_ds3_full_tail.ps1','notes/stage4-ds3-full-tail.md','literature/stage4-ds3-full-tail-sources.md',
            'manuscript/ds3-full-tail/main.tex','output/pdf/de-sitter-ds3-full-tail.pdf','sources/manifest.json']
    paths += [p.relative_to(root).as_posix() for p in sorted(base.iterdir()) if p.is_file()]
    evidence=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),status='stage verified; project integration pending; unfrozen',
                  pages=4,independent_checks=407,compilers=['two local pdfLaTeX passes','desktop compiler success'],
                  visual_review=dict(pages=[1,2,3,4],longest_edge_pixels=1100,status='all inspected without layout defects'),
                  pending_previous_stage=dict(session_id=47841,output='results/stage4/ds3-collision/controlled-validation.json',counted_as_passed=False),
                  hashes={p:sha(root/p) for p in paths})
    target.write_text(json.dumps(evidence,indent=2),encoding='utf-8')
    for p,digest in evidence['hashes'].items():assert sha(root/p)==digest
    print(json.dumps(dict(status='passed',bound_files=len(paths),old_files_unchanged=len(unchanged),old_zip_count=len(archives),changed_entries=changed),indent=2))


if __name__=='__main__':main()
