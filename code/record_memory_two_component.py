"""Bind the two-component memory proof, tests, diagnostics and final PDF."""
from pathlib import Path
from decimal import Decimal
import datetime
import hashlib
import json
import platform
import re
import numpy
import mpmath

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-two-component'


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


if __name__=='__main__':
    dest=OUT/'local-evidence.json'
    if dest.exists():raise FileExistsError(dest)
    mutable={'README.md','PLAN.md','NEXT.md'}
    parents={}
    for stage in ['susskind-memory-uniform-bound','susskind-memory-return-background',
                  'susskind-memory-charged-transfer','susskind-memory-wick-rank',
                  'susskind-memory-central-spectrum','susskind-memory-central-raw',
                  'susskind-memory-mixed-spectrum','susskind-memory-charge-sum']:
        name=f'results/{stage}/local-evidence.json'
        path=ROOT/name;record=json.loads(path.read_text(encoding='utf-8'));checked=[]
        for item,digest in record['hashes'].items():
            if item in mutable:continue
            assert sha(ROOT/item)==digest,(stage,item)
            checked.append(item)
        parents[name]=dict(sha256=sha(path),unchanged_scientific_files=checked)
    data={name:json.loads((OUT/name).read_text(encoding='utf-8')) for name in
          ['four-copy-q-checks.json','four-copy-finite-time-checks.json','neutral-coefficient-checks.json']}
    verified={}
    for name,record in data.items():
        for item,digest in record['input_sha256'].items():
            assert sha(ROOT/item)==digest,(name,item)
        verified[name]=record['input_sha256']
    qrows=data['four-copy-q-checks.json']['rows']
    totals={key:sum(row[key] for row in qrows) for key in
            ['conditional_trace_checks','three_color_replica_checks','OU_Wick_degrees','Taylor_checks']}
    assert len(qrows)==6 and list(totals.values())==[2520,7200,66,42]
    integrals=data['four-copy-finite-time-checks.json']['records']
    assert len(integrals)==24
    assert all(row['time_order_difference']<2e-10
               and row['integral_operator_norm']<=row['analytic_operator_bound']+2e-12
               and row['min_spectral_gap']>1e-8 for row in integrals)
    assert max(row['pinched_residual_HS_norm'] for row in integrals)>1e-5
    neutral=data['neutral-coefficient-checks.json']
    assert len(neutral['finite_CAR'])==12 and len(neutral['rank_slopes'])==360
    max_car=max(max(row[key] for key in ['first_conserved_difference_error',
                    'second_conserved_difference_error','exact_rank_covariance_error'])
                for row in neutral['finite_CAR'])
    assert max_car<3e-12
    assert all(Decimal(row['precision_difference'])<Decimal('1e-70')
               for row in neutral['rank_slopes'])
    build=ROOT/'tmp/pdfs/susskind-memory-two-component-build'
    pdf=ROOT/'output/pdf/de-sitter-susskind-memory-two-component.pdf'
    log=(OUT/'main.log').read_text(errors='replace')
    assert re.search(r'\(8 pages,',log)
    assert not re.search(r'undefined|Overfull|Rerun to get|Fatal error',log)
    assert sha(pdf)==sha(build/'main.pdf')
    assert sha(OUT/'main.log')==sha(build/'main.log')
    pngs=[build/f'page-{i}.png' for i in range(1,9)]
    passes=[build/f'pass{i}.txt' for i in [1,2]]
    assert all(path.exists() for path in pngs+passes)
    assert all(path.stat().st_mtime>=pdf.stat().st_mtime for path in pngs)
    names=['code/memory_four_copy_background.py','code/check_four_copy_finite_time.py',
           'code/check_neutral_memory_coefficients.py','code/build_memory_two_component.ps1',
           'code/record_memory_two_component.py',
           'manuscript/susskind-memory-two-component/main.tex',
           'output/pdf/de-sitter-susskind-memory-two-component.pdf',
           'notes/susskind-memory-two-component.md','notes/susskind-memory-neutral-series-seed.md',
           'literature/susskind-memory-two-component-sources.md',
           'journal/2026-11-07-susskind-memory-two-component.md',
           'README.md','PLAN.md','NEXT.md']
    paths=[ROOT/name for name in names]+[OUT/name for name in data]+[OUT/'main.log']
    paths+=sorted((OUT/'diagnostics').glob('*'))
    result=dict(
        status='Local research progress; full objective remains active',
        created_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        scientific_status='Four-u conditional force term has a strong leading energy regression and vanishes after collective energy pinching. Actual memory field reduces to charged defect A plus neutral residual B in full-trace L2 on compact real memory times; the covariance is their full four-term sum. First two neutral Taylor coefficients and nonzero second-order cross measures identified. The joint finite-time kernel and bulk dictionary remain open.',
        git_commit_created=False,git_push_performed=False,archive_created=False,external_reply_sent=False,
        previous_evidence=parents,verified_scientific_inputs=verified,
        checks=dict(q_parameter_sets=6,exact_q_equalities=totals,finite_CAR_integrals=24,
                    max_integral_order_difference=max(x['time_order_difference'] for x in integrals),
                    max_finite_N_pinched_norm=max(x['pinched_residual_HS_norm'] for x in integrals),
                    CAR_coefficient_cases=12,max_CAR_coefficient_error=max_car,rank_slope_cases=360,
                    max_relative_slope_error_at_p128=max(abs(float(x['relative_error'])) for x in neutral['rank_slopes'] if x['p']==128 and x['deficit'])),
        parameters='Analytic result: standard CAR, same original J across charges, overlapping complex p-body, full normalized beta=mu=0 trace, hbar=1/string real time, even N, fixed Js and compact memory times. Explicit spectra assume lambda_N->lambda>0; reduction also holds for positive compact lambda sequences. All finite fixtures and relative errors are in the scientific JSON.',
        environment=dict(python=platform.python_version(),numpy=numpy.__version__,mpmath=mpmath.__version__),
        evidence_scopes={name:record['scope'] for name,record in data.items()},
        proof_scope='Fixed colored chord moments and the q-Fock commuting square yield strong polynomial regression; a total-degree Gaussian moment tail gives uniform real time. Bounded extension, alpha drift and binomial tails imply weighted pinching. Exact neutral coefficient identities and existing marked spectra supply the cross calibration. No derivative exchange with an unproved finite-time limit is claimed.',
        literature_scope='BKS conditional expectations and q-Wick/Hermite structure are established results. Colored extension and the specified finite-memory reduction are derived applications; no originality priority claim.',
        continuation_seed_status='Neutral-series exact algebraic rearrangement and scaled Jacobi error recurrence saved as next calculation interface. No new all-degree joint fluctuation limit is claimed there; independent finite-polynomial check and rank-tail control remain next.',
        build=dict(command='powershell -NoProfile -ExecutionPolicy Bypass -File code/build_memory_two_component.ps1',
                   pages=8,passes=2,compiler_header=log.splitlines()[0],
                   render_command='pdftoppm -r 95 -png output/pdf/de-sitter-susskind-memory-two-component.pdf tmp/pdfs/susskind-memory-two-component-build/page',
                   final_pages_visually_inspected=list(range(1,9)),
                   visual_result='Final eight pages inspected after the last build; readable equations, correct references, clean margins and no clipping or overlaps.'),
        diagnostics=[
            'The first finite-time driver failed its final nonconservation-witness assertion after all integral tests passed; no JSON was written. Initial p=1/edge sectors were too symmetric. Source and failure text retained; adding the interior N5p2m3 fixture supplied the nonzero witness.',
            'The q exact driver and neutral coefficient driver passed first execution. Existing scientific outputs were not overwritten.',
            'First LaTeX build failed because AA is a built-in macro. Second build had two overfull boxes. Sources and logs retained; renaming the macro and wrapping equations fixed them.',
            'Read-only mathematical review added positive-part kappa_2p definition, the linear-Gaussian multiplication bound and the scope of pinching for energy descendants. Final two-pass compile and all-page rendering passed.'],
        process_state='All scientific, build and rendering commands terminal; no live calculation handles.',
        hashes={path.relative_to(ROOT).as_posix():sha(path) for path in paths},
        render_hashes={path.relative_to(ROOT).as_posix():sha(path) for path in pngs},
        build_pass_hashes={path.relative_to(ROOT).as_posix():sha(path) for path in passes})
    with dest.open('x',encoding='utf-8') as stream:stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(dict(status='passed',pages=8,stage_files=len(paths),
        previous_scientific_files_checked=sum(len(x['unchanged_scientific_files']) for x in parents.values()),
        exact_q_checks=sum(totals.values()),finite_integrals=24,coefficient_cases=12,slope_cases=360)))
