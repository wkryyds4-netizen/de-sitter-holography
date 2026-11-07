"""Read-only validation of scientific evidence, then save a local manifest."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import platform
import re
import subprocess

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / 'results/susskind-memory-raw-spectrum'
BUILD = ROOT / 'tmp/pdfs/susskind-memory-raw-spectrum-build'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    target = OUT / 'local-evidence.json'
    if target.exists():
        raise FileExistsError(target)
    records = {}
    for stem in ['symbolic-spectrum', 'finite-restoration', 'kernel-checks']:
        data = json.loads((OUT / (stem + '.json')).read_text())
        assert data['status'] == 'passed'
        for name, digest in data['input_sha256'].items():
            assert sha(ROOT / name) == digest, name
        records[stem] = data
    finite = records['finite-restoration']['records']
    assert len(finite) == 4
    assert sum(r['exact_identities'] for r in finite) == 16
    assert len(records['symbolic-spectrum']['unit_variance_probe_moments']) == 8
    kernels = records['kernel-checks']['kernel']['records']
    assert [r['cutoff'] for r in kernels] == [38, 57, 146, 314]
    assert all(r['normalized_raw_truncation_bound'] < 1e-12 for r in kernels)
    parentpath = ROOT / 'results/susskind-memory-mixed-spectrum/local-evidence.json'
    parent = json.loads(parentpath.read_text())
    unchanged = []
    for name, digest in parent['hashes'].items():
        if name in ['README.md', 'PLAN.md', 'NEXT.md']:
            continue
        assert sha(ROOT / name) == digest, name
        unchanged.append(name)
    log = (OUT / 'main.log').read_text(errors='replace')
    assert re.search(r'6\s+pages', log)
    assert not re.search(r'Overfull|undefined|Rerun to get|Fatal error', log)
    tex = ROOT / 'manuscript/susskind-memory-raw-spectrum/main.tex'
    assert not any(ord(c) < 32 and c not in '\n\r' for c in tex.read_text())
    pdf = ROOT / 'output/pdf/de-sitter-susskind-memory-raw-spectrum.pdf'
    assert sha(pdf) == sha(BUILD / 'main.pdf')
    info = subprocess.check_output(['pdfinfo', str(pdf)], text=True, errors='replace')
    assert re.search(r'Pages:\s+6\b', info)
    renders = [BUILD / f'page-{i}.png' for i in range(1, 7)]
    assert all(p.exists() for p in renders)
    head = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
    assert head == parent['current_git_head']
    names = ['code/memory_raw_spectrum.py', 'code/check_memory_raw_restoration.py',
             'code/memory_raw_kernel.py', 'code/build_memory_raw_spectrum.ps1',
             'code/record_memory_raw_spectrum.py',
             'manuscript/susskind-memory-raw-spectrum/main.tex',
             'output/pdf/de-sitter-susskind-memory-raw-spectrum.pdf',
             'notes/susskind-memory-raw-spectrum.md',
             'literature/susskind-memory-raw-spectrum-sources.md',
             'journal/2026-11-06-susskind-memory-raw-spectrum.md',
             'README.md', 'PLAN.md', 'NEXT.md']
    paths = [ROOT / n for n in names] + sorted(p for p in OUT.rglob('*') if p.is_file())
    data = dict(
        created_at_utc=datetime.now(timezone.utc).isoformat(),
        status='Saved locally; no publication',
        scientific_status='At fixed half filling, restore dynamic descendants of the raw second memory coefficient after an explicit conserved subtraction. Complete two-energy self and cross spectra assembled; ordered small-lambda normalized shape is (9 sech^8 - 2 sech^6)/7. Uniform central-charge extension, finite-memory remainder and bulk dictionary remain open.',
        current_git_head=head, git_commit_created=False, git_push_performed=False, archive_created=False,
        parent_manifest_sha256=sha(parentpath), previous_stage_scientific_files_unchanged=unchanged,
        checks=dict(exact_finite_identities=16, finite_systems=4, probe_moments=8,
                    kernel_cutoffs=[38, 57, 146, 314], refined_levels_added=12),
        environment=dict(python=platform.python_version(), calculations=[r['environment'] for r in records.values()]),
        build=dict(command='powershell -NoProfile -ExecutionPolicy Bypass -File code/build_memory_raw_spectrum.ps1',
                   compiler=log.splitlines()[0], passes_per_successful_build=2, successful_builds=1, pages=6,
                   render_command='pdftoppm -r 95 -png output/pdf/de-sitter-susskind-memory-raw-spectrum.pdf tmp/pdfs/susskind-memory-raw-spectrum-build/page',
                   visual_review='All six final page PNGs viewed at 95 dpi; equations, citations, table and page boundaries clear.'),
        failures=[],
        process_state='Scientific calculations, compilation and rendering terminal; no active calculation.',
        evidence_scope='Exact algebra checks and weighted spectral checks supplement established asymptotic input proofs. Finite CAR checks do not establish finite-N spectral equality; floating checks are not interval certification. Second Taylor polynomial is not a finite-memory kernel.',
        hashes={p.relative_to(ROOT).as_posix(): sha(p) for p in paths},
        render_hashes={p.relative_to(ROOT).as_posix(): sha(p) for p in renders},
        build_pass_hashes={f'pass{i}.txt': sha(BUILD / f'pass{i}.txt') for i in [1, 2]})
    target.write_text(json.dumps(data, indent=2) + '\n')
    print(json.dumps(dict(status='saved locally', pages=6, files=len(paths),
                         prior_scientific_files_unchanged=len(unchanged), git_head=head)))


if __name__ == '__main__':
    main()
