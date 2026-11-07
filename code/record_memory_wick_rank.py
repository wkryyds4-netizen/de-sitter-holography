"""Bind the all-chaos rank stage without rerunning or replacing scientific outputs."""
from pathlib import Path
import datetime
from decimal import Decimal
import hashlib
import json
import re

ROOT = Path(__file__).resolve().parent.parent
STAGE = ROOT / 'results/susskind-memory-wick-rank'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


if __name__ == '__main__':
    dest = STAGE / 'local-evidence.json'
    if dest.exists():
        raise FileExistsError(dest)
    parent = ROOT / 'results/susskind-memory-return-background/local-evidence.json'
    previous = json.loads(parent.read_text(encoding='utf-8'))
    checked_previous = []
    for name, digest in previous['hashes'].items():
        if name in ('README.md', 'PLAN.md', 'NEXT.md'):
            continue
        assert sha(ROOT / name) == digest, name
        checked_previous.append(name)

    names = ('rank-chain-checks.json', 'chaos-and-probe-checks.json',
             'creation-commutator-checks.json')
    outputs = {name: json.loads((STAGE / name).read_text(encoding='utf-8'))
               for name in names}
    input_verification = {}
    for name, data in outputs.items():
        for source, digest in data['input_sha256'].items():
            assert sha(ROOT / source) == digest, (name, source)
        input_verification[name] = data['input_sha256']
    rank, interface, creation = (outputs[name] for name in names)
    assert rank['exact']['exact_equalities'] == 252
    assert rank['rational_bound_instances'] == 1800
    assert rank['symbolic']['general_rational_identities'] == 4
    assert len(rank['symbolic']['fixed_k_symbolic_checks']) == 6
    assert len(rank['CAR_channels']) == 3
    assert len(rank['high_precision']) == 5
    assert all(Decimal(row['precision_difference']) < Decimal('1e-60')
               for row in rank['high_precision'])
    assert len(interface['highest_chaos']) == 9
    assert len(interface['charged_probes']) == 3
    assert sum(len(system['cases']) for system in creation['systems']) == 18
    assert creation['exact_step_cases'] == 43

    render = ROOT / 'tmp/pdfs/susskind-memory-wick-rank-build'
    log = (STAGE / 'main.log').read_text(errors='replace')
    assert not re.search(r'Overfull|undefined|Rerun to get|Fatal error', log)
    assert re.search(r'\(8 pages,', log)
    pdf = ROOT / 'output/pdf/de-sitter-susskind-memory-wick-rank.pdf'
    assert sha(pdf) == sha(render / 'main.pdf')
    assert sha(STAGE / 'main.log') == sha(render / 'main.log')
    pngs = [render / f'page-{i}.png' for i in range(1, 9)]
    passes = [render / f'pass{i}.txt' for i in (1, 2)]
    assert all(path.exists() for path in pngs + passes)
    assert all(path.stat().st_mtime >= pdf.stat().st_mtime for path in pngs)

    paths = [ROOT / name for name in (
        'code/memory_wick_rank.py', 'code/check_memory_wick_interface.py',
        'code/check_memory_chaos_creation.py', 'code/build_memory_wick_rank.ps1',
        'code/record_memory_wick_rank.py',
        'manuscript/susskind-memory-wick-rank/main.tex',
        'output/pdf/de-sitter-susskind-memory-wick-rank.pdf',
        'notes/susskind-memory-wick-rank.md',
        'literature/susskind-memory-wick-rank-sources.md',
        'journal/2026-11-07-susskind-memory-wick-rank.md',
        'README.md', 'PLAN.md', 'NEXT.md')]
    paths += [STAGE / name for name in names] + [STAGE / 'main.log']
    paths += sorted((STAGE / 'diagnostics').glob('*'))
    result = dict(
        status='Saved locally; overall research goal remains active',
        created_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        scientific_status=(
            'Uniform rank-chaos deficit bound for the full ordinary exp(isH) on '
            'compact real-time intervals, including lower chaoses. The connection '
            'to actual adjacent-charge force and energy regression, X0 tightness, '
            'finite-memory covariance and the bulk dictionary remains open.'),
        git_commit_created=False, git_push_performed=False, archive_created=False,
        external_reply_sent=False,
        previous_evidence={parent.relative_to(ROOT).as_posix(): dict(
            sha256=sha(parent), unchanged_scientific_files=checked_previous)},
        verified_scientific_inputs=input_verification,
        checks=dict(exact_rank_equalities=252, rational_bound_instances=1800,
                    general_symbolic_identities=4, fixed_k_symbolic_checks=6,
                    independent_CAR_channels=3, high_precision_rows=5,
                    highest_chaos_tensors=9, charged_probe_systems=3,
                    full_Hermite_cases=18, exact_one_step_laws=43),
        environments={name: data['environment'] for name, data in outputs.items()},
        parameters_location='Each scientific JSON contains its finite systems and parameters; charged-probe seed is 1107001.',
        numerical_scopes={name: data['scope'] for name, data in outputs.items()},
        proof_scope=(
            'The uniform theorem is the analytic argument in manuscript Sections '
            '2-6, with two independent read-only collaborator checks. Finite '
            'coefficient and rational checks are not substitutes for this proof. '
            'The rank driver predates the extension to the ordinary propagator; '
            'its narrower saved scope is retained unchanged. Large k=100 tests '
            'Casimir evaluation, not a fixed-k limit; precision comparisons are '
            'not interval-rounding certificates. Real time only.'),
        build=dict(
            command='powershell -NoProfile -ExecutionPolicy Bypass -File code/build_memory_wick_rank.ps1',
            passes=2, pages=8, compiler_log_header=log.splitlines()[0],
            render_command='pdftoppm -r 95 -png output/pdf/de-sitter-susskind-memory-wick-rank.pdf tmp/pdfs/susskind-memory-wick-rank-build/page',
            final_pages_visually_inspected=list(range(1, 9)),
            visual_result='No clipped formulas, overlaps, broken references or literal qquad remains.'),
        diagnostics=[
            'First six-page draft compiled cleanly but visual review found literal qquad. Its source and PDF are preserved in diagnostics; final source corrected and expanded to eight pages.',
            'Old PDF skill 26.1030.10813 marker path failed MODULE_NOT_FOUND after runtime update. Skill 26.1105.12137 was read and marker succeeded exactly once before authoring.',
            'All three scientific drivers passed on their first run. Existing output files were not regenerated.',
            'Read-only collaborator probes encountered series-truncation and console-encoding errors, corrected without changing scientific evidence.'
        ],
        process_state='All scientific and build processes terminal; no active shell handle.',
        hashes={path.relative_to(ROOT).as_posix(): sha(path) for path in paths},
        render_hashes={path.relative_to(ROOT).as_posix(): sha(path) for path in pngs},
        build_pass_hashes={path.relative_to(ROOT).as_posix(): sha(path) for path in passes})
    with dest.open('x', encoding='utf-8') as stream:
        stream.write(json.dumps(result, indent=2, ensure_ascii=False) + '\n')
    print(json.dumps(dict(status='passed', pages=8,
                         previous_files_checked=len(checked_previous),
                         scientific_outputs_checked=len(outputs),
                         saved_files=len(paths))))
