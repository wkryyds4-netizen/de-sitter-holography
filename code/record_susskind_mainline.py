"""Bind the reviewed mainline programme; no spectral computation is claimed."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import re
import subprocess

ROOT = Path(__file__).resolve().parent.parent
DEST = ROOT / 'results/susskind-mainline'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def command_output(args):
    result = subprocess.run(args, capture_output=True, check=True)
    return (result.stdout + result.stderr).decode('utf-8', errors='replace')


def save(name, data):
    with (DEST / name).open('x', encoding='utf-8') as handle:
        json.dump(data, handle, ensure_ascii=False, indent=2)
        handle.write('\n')


def main():
    for name in ['build-review.json', 'source-record.json', 'evidence.json']:
        if (DEST / name).exists():
            raise FileExistsError(DEST / name)
    log = (DEST / 'main.log').read_text(encoding='utf-8', errors='replace')
    assert not re.search(r'undefined|Rerun to get|Overfull|Fatal error', log)
    assert '(7 pages,' in log
    info = command_output(['pdfinfo', str(ROOT / 'output/pdf/de-sitter-susskind-mainline.pdf')])
    assert re.search(r'Pages:\s+7\b', info)
    rendered = [ROOT / f'tmp/pdfs/susskind-mainline-build/page-{i}.png' for i in range(1, 8)]
    assert all(p.is_file() for p in rendered)
    now = datetime.now(timezone.utc).isoformat()
    save('build-review.json', {
        'created_at_utc': now, 'pages': 7, 'latex_passes': 2,
        'compiler': command_output(['pdflatex', '--version']).splitlines()[0],
        'renderer': command_output(['pdftoppm', '-v']).splitlines()[0],
        'pdfinfo': info, 'overflow_or_undefined_references': False,
        'visual_review': 'Assistant inspected all seven final 1200-pixel PNG pages: no clipping, overlap or broken formula layout.',
        'rendered_pages_local': {p.relative_to(ROOT).as_posix(): digest(p) for p in rendered},
        'build_command': 'powershell -File code/build_susskind_mainline.ps1',
    })
    ledger = ROOT / 'literature/susskind-mainline.md'
    save('source-record.json', {
        'created_at_utc': now, 'ledger': ledger.relative_to(ROOT).as_posix(),
        'urls': sorted(set(re.findall(r'https://arxiv\.org/(?:abs|html)/[^)\s]+', ledger.read_text(encoding='utf-8')))),
        'reading_scope': 'Passages and unread follow-ups are recorded individually in the ledger. Online readings, not a claim to have reproduced every cited derivation.',
        'new_fulltext_cache': None,
    })
    paths = [
        'code/build_susskind_mainline.ps1', 'code/record_susskind_mainline.py',
        'manuscript/susskind-mainline/main.tex', 'output/pdf/de-sitter-susskind-mainline.pdf',
        'notes/susskind-mainline.md', 'literature/susskind-mainline.md',
        'results/susskind-mainline/main.log', 'results/susskind-mainline/build-review.json',
        'results/susskind-mainline/source-record.json',
    ]
    save('evidence.json', {
        'created_at_utc': now, 'status': 'programme and finite-system identities recorded',
        'scientific_status': 'Literature-led construction; no new numerical spectrum or established bulk dual.',
        'hashes': {name: digest(ROOT / name) for name in paths},
    })
    print(json.dumps({'files_bound': len(paths), 'pages': 7, 'status': 'saved'}))


if __name__ == '__main__':
    main()
