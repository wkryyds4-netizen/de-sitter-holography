"""Bind the local conserved-return stage and verify previous evidence bytes."""
from pathlib import Path
import datetime
import hashlib
import json
import re

ROOT = Path(__file__).resolve().parent.parent
STAGE = ROOT / 'results/susskind-memory-return-background'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


if __name__ == '__main__':
    dest = STAGE / 'local-evidence.json'
    if dest.exists():
        raise FileExistsError(dest)
    parents = [ROOT / 'results/susskind-memory-return-expansion/local-evidence.json',
               ROOT / 'results/issue-1-review/local-evidence.json']
    previous = {}
    for parent in parents:
        record = json.loads(parent.read_text(encoding='utf-8'))
        names = []
        for name, digest in record['hashes'].items():
            if name in ['README.md', 'PLAN.md', 'NEXT.md']:
                continue
            assert sha(ROOT / name) == digest, name
            names.append(name)
        previous[parent.relative_to(ROOT).as_posix()] = dict(
            sha256=sha(parent), unchanged_scientific_files=names)
    check = json.loads((STAGE / 'bubble-and-return-checks.json').read_text())
    log = (STAGE / 'main.log').read_text(errors='replace')
    assert not re.search(r'Overfull|undefined|Rerun to get|Fatal error', log)
    assert re.search(r'\(6 pages,', log)
    assert check['chords']['total_exact_equalities'] == 2373
    assert len(check['return_background']) == 27
    paths = [ROOT / x for x in [
        'code/memory_return_background.py', 'code/build_memory_return_background.ps1',
        'code/record_memory_return_background.py',
        'manuscript/susskind-memory-return-background/main.tex',
        'output/pdf/de-sitter-susskind-memory-return-background.pdf',
        'notes/susskind-memory-return-background.md',
        'literature/susskind-memory-return-background-sources.md',
        'journal/2026-11-06-susskind-memory-return-background.md',
        'README.md', 'PLAN.md', 'NEXT.md']]
    paths += list(STAGE.glob('*.json')) + list(STAGE.glob('*.log'))
    paths += list((STAGE / 'diagnostics').glob('*'))
    render = ROOT / 'tmp/pdfs/susskind-memory-return-background-build'
    pngs = [render / f'page-{i}.png' for i in range(1, 7)]
    assert all(p.exists() for p in pngs)
    result = dict(
        status='Saved locally; overall research goal remains active',
        created_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        scientific_status='Strong leading ordinary-kernel approximation and conserved first-return limit. Pinched field L2 equivalence; uniform L2 boundedness and finite-memory covariance remain open.',
        git_commit_created=False, git_push_performed=False, archive_created=False,
        external_reply_sent=False, previous_evidence=previous,
        checks=dict(exact_fixed_bubble_equalities=2373, scalar_values=27,
                    max_scalar_order_difference=max(x['absolute_order_difference'] for x in check['return_background'])),
        environment=check['environment'], parameters=check['parameters'],
        build=dict(command='powershell -NoProfile -ExecutionPolicy Bypass -File code/build_memory_return_background.ps1',
                   pages=6, passes=2,
                   render_command='pdftoppm -r 95 -png output/pdf/de-sitter-susskind-memory-return-background.pdf tmp/pdfs/susskind-memory-return-background-build/page',
                   all_six_final_pages_visually_inspected=True),
        failures=['First layout had three overfull lines; diagnostics/first-layout.log preserved. Scientific driver passed its first run.',
                  'BKS web find requests failed after a successful open; earlier source reading and collaborator theorem-location checks retained.',
                  'Initial literature-note Xu title mismatch corrected against arXiv abstract before final save.'],
        evidence_scope=check['scope'],
        process_state='All current scientific and build processes terminal; no active shell handle.',
        hashes={p.relative_to(ROOT).as_posix(): sha(p) for p in paths},
        render_hashes={p.relative_to(ROOT).as_posix(): sha(p) for p in pngs},
        build_pass_hashes={p.relative_to(ROOT).as_posix(): sha(p) for p in render.glob('pass*.txt')})
    dest.write_text(json.dumps(result, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    print(json.dumps(dict(status='passed', pages=6,
                         previous_files_checked=sum(len(x['unchanged_scientific_files']) for x in previous.values()))))
