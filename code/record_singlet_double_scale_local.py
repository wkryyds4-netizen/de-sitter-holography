"""Bind the completed local stage without committing or publishing it."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import re
import subprocess

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / 'results/susskind-singlet-double-scale'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    target = OUT / 'local-evidence.json'
    if target.exists():
        raise FileExistsError(target)
    previous_path = ROOT / 'results/susskind-singlet-moments/local-evidence.json'
    previous = json.loads(previous_path.read_text(encoding='utf-8'))
    for path, digest in previous['hashes'].items():
        assert sha(ROOT / path) == digest, path
    for name in ['validation.json', 'centered-gram.json']:
        data = json.loads((OUT / name).read_text(encoding='utf-8'))
        assert data['status'].startswith('passed'), name
        for path, digest in data['input_sha256'].items():
            assert sha(ROOT / path) == digest, path
    log = (OUT / 'main.log').read_text(encoding='utf-8', errors='replace')
    assert not re.search(r'Overfull|undefined|Fatal error|Rerun to get', log)
    assert '(7 pages,' in log
    build_dir = ROOT / 'tmp/pdfs/susskind-singlet-double-scale-build'
    for i in [1, 2]:
        raw = (build_dir / f'pass{i}.txt').read_bytes()
        encoding = 'utf-16' if raw.startswith((b'\xff\xfe', b'\xfe\xff')) else 'utf-8'
        assert '(7 pages,' in raw.decode(encoding, errors='replace')
    assert sha(build_dir / 'main.pdf') == sha(ROOT / 'output/pdf/de-sitter-susskind-singlet-double-scale.pdf')
    renders = [build_dir / f'page-{i}.png' for i in range(1, 8)]
    assert all(p.exists() for p in renders)
    files = [ROOT / p for p in [
        'code/singlet_coefficient.py', 'code/check_singlet_coefficient.py',
        'code/run_singlet_double_scale.py', 'code/run_singlet_centered_moment.py',
        'code/finish_singlet_double_scale.py', 'code/build_singlet_double_scale.ps1',
        'code/record_singlet_double_scale_local.py',
        'manuscript/susskind-singlet-double-scale/main.tex',
        'output/pdf/de-sitter-susskind-singlet-double-scale.pdf',
        'notes/susskind-singlet-double-scale.md',
        'literature/susskind-singlet-double-scale-sources.md',
        'journal/2026-11-05-susskind-singlet-double-scale.md',
        'README.md', 'PLAN.md', 'NEXT.md',
    ]]
    files += sorted(p for p in OUT.iterdir() if p.is_file())
    payload = dict(
        created_at_utc=datetime.now(timezone.utc).isoformat(),
        status='local increment recorded; not frozen or published',
        scientific_status='Double-scaled equal-time variance and second frequency moment; full spectrum and bulk dictionary remain open.',
        current_git_head=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        git_commit_created=False, git_push_performed=False, archive_created=False,
        checks='Prior-stage hashes and current validation input hashes verified; final build checked. No global project audit or archive rehash.',
        previous_stage=dict(record_sha256=sha(previous_path), bound_files_verified=len(previous['hashes'])),
        build=dict(
            command='powershell -NoProfile -ExecutionPolicy Bypass -File code/build_singlet_double_scale.ps1',
            compiler=subprocess.check_output(['pdflatex', '--version'], text=True, errors='replace').splitlines()[0],
            passes=2, pages=7,
            visual_review='All seven 1200-pixel PNG pages inspected earlier in this research turn; no clipping or overlap.',
            initial_issue='Long filename overflow repaired; initial log retained.',
        ),
        hashes={p.relative_to(ROOT).as_posix(): sha(p) for p in files},
        local_render_hashes={p.relative_to(ROOT).as_posix(): sha(p) for p in renders},
        build_pass_hashes={f'pass{i}.txt': sha(build_dir / f'pass{i}.txt') for i in [1, 2]},
    )
    with target.open('x', encoding='utf-8') as f:
        json.dump(payload, f, indent=2)
        f.write('\n')
    print(json.dumps(dict(status='recorded locally', files=len(files), previous_files_verified=len(previous['hashes']), pages=7, git_commit=False, git_push=False)))


if __name__ == '__main__':
    main()
