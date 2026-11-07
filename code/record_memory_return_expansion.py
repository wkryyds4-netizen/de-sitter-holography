"""Bind the locally completed return-expansion stage; never publish or overwrite."""
from pathlib import Path
import json,hashlib,datetime,re

ROOT=Path(__file__).resolve().parent.parent
STAGE=ROOT/'results/susskind-memory-return-expansion'
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()

if __name__=='__main__':
    dest=STAGE/'local-evidence.json'
    if dest.exists(): raise FileExistsError(dest)
    parent=ROOT/'results/susskind-memory-charge-sum/local-evidence.json'
    previous=json.loads(parent.read_text(encoding='utf-8'))
    unchanged=[]
    for name,digest in previous['hashes'].items():
        if name in ['README.md','PLAN.md','NEXT.md']: continue
        assert sha(ROOT/name)==digest,name
        unchanged.append(name)
    check=json.loads((STAGE/'finite-checks.json').read_text())
    log=(STAGE/'main.log').read_text(errors='replace')
    assert not re.search(r'Overfull|undefined|Rerun to get|Fatal error',log)
    assert re.search(r'\(6 pages,',log)
    paths=[ROOT/x for x in [
        'code/memory_return_expansion.py','code/build_memory_return_expansion.ps1',
        'code/record_memory_return_expansion.py','manuscript/susskind-memory-return-expansion/main.tex',
        'output/pdf/de-sitter-susskind-memory-return-expansion.pdf',
        'notes/susskind-memory-return-expansion.md',
        'literature/susskind-memory-return-expansion-sources.md',
        'journal/2026-11-06-susskind-memory-return-expansion.md']]
    paths+=list(STAGE.glob('*.json'))+list(STAGE.glob('*.log'))+list((STAGE/'diagnostics').glob('*'))
    render=ROOT/'tmp/pdfs/susskind-memory-return-expansion-build'
    pngs=[render/f'page-{i}.png' for i in range(1,7)]
    assert all(p.exists() for p in pngs)
    result=dict(status='Saved locally; goal remains active',
        created_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        scientific_status='Finite-return field O(p^-2) and absolute covariance O(p^-1) reduction. Finite-memory tightness and limiting kernel still to compute.',
        git_commit_created=False,git_push_performed=False,archive_created=False,
        parent_manifest_sha256=sha(parent),previous_scientific_files_unchanged=unchanged,
        checks=dict(force_diagonal_entries=sum(x['diagonal_entries_checked'] for x in check['force_square_functions']['exact_CAR']),
                    rational_bounds=check['force_square_functions']['rational_bound_checks'],
                    finite_remainder_bounds=48,common_coupling_CAR_systems=4,scalar_cases=5),
        environment=check['environment'],
        build=dict(command='powershell -NoProfile -ExecutionPolicy Bypass -File code/build_memory_return_expansion.ps1',
                   pages=6,passes=2,render_command='pdftoppm -r 95 -png output/pdf/de-sitter-susskind-memory-return-expansion.pdf tmp/pdfs/susskind-memory-return-expansion-build/page',
                   all_six_pages_visually_inspected=True),
        failures=['First layout had three overfull lines; first-layout.log preserved. Scientific driver passed on its first run.'],
        evidence_scope='Exact CAR covariance counts and finite-matrix identities; asymptotic rates follow the written inequalities, not finite-size fits.',
        process_state='Calculation and final build terminal; no active return-expansion process.',
        hashes={p.relative_to(ROOT).as_posix():sha(p) for p in paths},
        render_hashes={p.relative_to(ROOT).as_posix():sha(p) for p in pngs},
        build_pass_hashes={p.relative_to(ROOT).as_posix():sha(p) for p in render.glob('pass*.txt')})
    dest.write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(json.dumps(dict(status='passed',previous_files_unchanged=len(unchanged),pages=6)))
