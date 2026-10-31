"""Bind the peak-stage record without overwriting previous scientific evidence."""
import hashlib
import json
import zipfile
from datetime import datetime, timezone
from pathlib import Path


def main():
    root = Path(__file__).resolve().parents[1]
    base = root / 'results/stage4/ds3-peak-limit'
    dest = base / 'evidence.json'
    if dest.exists() or (base / 'preservation.json').exists():
        raise RuntimeError('Preserve existing evidence')
    sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    prior = {}
    for stage in ['ds3-full-tail', 'ds3-shadow-duality']:
        old = json.loads((root / f'results/stage4/{stage}/evidence.json').read_text())
        for name, digest in old['hashes'].items():
            assert sha(root / name) == digest, name
        prior[stage] = len(old['hashes'])
    study = json.loads((base / 'study.json').read_text())
    audit = json.loads((base / 'audit-refined.json').read_text())
    for record in [study, audit]:
        for name, digest in record['hashes'].items():
            assert sha(root / name) == digest, name
    assert study['status'] == audit['status'] == 'passed'
    assert len(study['rows']) == 12 and audit['checks'] == 191
    assert json.loads((base / 'audit.json').read_text())['status'] == 'failed'
    baseline = root / 'snapshots/2026-10-31-ds3-bose.zip'
    assert sha(baseline) == 'fddbeeb767dca73e38b352c956be56d109e9f7283c436c51c100de31ba11f74c'
    allowed = {'README.md', 'PLAN.md', 'NEXT.md', 'journal/2026-10-31.md',
               'results/project-audit.json', 'results/project-files.json'}
    same, changed = 0, []
    with zipfile.ZipFile(baseline) as archive:
        for name in archive.namelist():
            if name.endswith('/'):
                continue
            if (root / name).is_file() and sha(root / name) == hashlib.sha256(archive.read(name)).hexdigest():
                same += 1
            else:
                assert name in allowed, name
                changed.append(name)
    archives = []
    for p in sorted((root / 'snapshots').glob('*.zip')):
        digest = sha(p)
        side = p.with_suffix('.sha256')
        if side.exists():
            assert side.read_text().split()[0] == digest, p
        archives.append(dict(file=p.name, sha256=digest, sidecar_present=side.exists()))
    assert len(archives) == 45
    preservation = dict(status='passed', old_snapshot_files_unchanged=same,
                        changed_entries=changed, prior_stage_bound_files_checked=prior,
                        old_archives=archives)
    (base / 'preservation.json').write_text(json.dumps(preservation, indent=2), encoding='utf-8')
    paths = ['code/' + n + '.py' for n in ['ds3_peak_normalization', 'run_ds3_peak_limit',
             'audit_ds3_peak_limit', 'audit_ds3_peak_limit_refined', 'record_ds3_peak_limit']]
    paths += ['code/build_ds3_peak_limit.ps1', 'notes/stage4-ds3-peak-limit.md',
              'literature/stage4-ds3-peak-limit-sources.md', 'manuscript/ds3-peak-limit/main.tex',
              'output/pdf/de-sitter-ds3-peak-limit.pdf', 'sources/manifest.json']
    paths += [p.relative_to(root).as_posix() for p in sorted(base.iterdir()) if p.is_file()]
    evidence = dict(created_at_utc=datetime.now(timezone.utc).isoformat(),
                    status='stage checks completed; project integration and freeze pending',
                    pages=5, study_runs=12, independent_audit_checks=191,
                    initial_failed_audit_preserved=True,
                    compilers=['local pdfLaTeX twice', 'desktop success'],
                    visual_review=dict(pages=[1, 2, 3, 4, 5], longest_edge_pixels=1100,
                                       status='all inspected, no layout defects'),
                    pending_prior_job=dict(session_id=47841, counted_as_passed=False),
                    hashes={p: sha(root / p) for p in paths})
    dest.write_text(json.dumps(evidence, indent=2), encoding='utf-8')
    for p, digest in evidence['hashes'].items():
        assert sha(root / p) == digest, p
    print(json.dumps(dict(status='passed', evidence_files=len(paths), old_files_unchanged=same,
                          prior_stage_bound_files=prior, zip_count=len(archives)), indent=2))


if __name__ == '__main__':
    main()
