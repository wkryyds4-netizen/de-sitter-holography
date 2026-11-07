"""Seal adjacent-charge propagation evidence, checking immutable predecessors."""
from pathlib import Path
import datetime,hashlib,json,re,platform
import numpy,scipy,sympy,mpmath

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-charged-transfer'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


if __name__=='__main__':
    dest=OUT/'local-evidence.json'
    if dest.exists():raise FileExistsError(dest)
    mutable={'README.md','PLAN.md','NEXT.md'}
    parents={}
    for name in ['results/susskind-memory-wick-rank/local-evidence.json',
                 'results/susskind-memory-wick-rank/continuation-seed.json']:
        path=ROOT/name;record=json.loads(path.read_text(encoding='utf-8'));checked=[]
        for item,digest in record['hashes'].items():
            if item in mutable:continue
            assert sha(ROOT/item)==digest,item
            checked.append(item)
        parents[name]=dict(sha256=sha(path),unchanged_scientific_files=checked)
    data={name:json.loads((OUT/name).read_text(encoding='utf-8')) for name in
          ['transfer-checks.json','polar-and-reverse-checks.json','closure-checks.json']}
    verified={}
    for name,result in data.items():
        for source,digest in result['input_sha256'].items():
            assert sha(ROOT/source)==digest,(name,source)
        verified[name]=result['input_sha256']
        assert len(result['systems'])==6
    assert data['transfer-checks.json']['symbolic']['exact_identities']==2
    assert data['polar-and-reverse-checks.json']['symbolic']['exact_squared_amplitude_identity']
    assert data['closure-checks.json']['exact_source_variance_cases']==230
    maxerrors={name:max(max(system['errors'].values()) for system in result['systems'])
               for name,result in data.items()}
    assert max(maxerrors.values())<3e-10
    build=ROOT/'tmp/pdfs/susskind-memory-charged-transfer-build'
    pdf=ROOT/'output/pdf/de-sitter-susskind-memory-charged-transfer.pdf'
    log=(OUT/'main.log').read_text(errors='replace')
    assert re.search(r'\(8 pages,',log)
    assert not re.search(r'undefined|Overfull|Rerun to get|Fatal error',log)
    assert sha(pdf)==sha(build/'main.pdf')
    assert sha(OUT/'main.log')==sha(build/'main.log')
    pngs=[build/f'page-{i}.png' for i in range(1,9)]
    passes=[build/f'pass{i}.txt' for i in (1,2)]
    assert all(x.exists() for x in pngs+passes)
    assert all(x.stat().st_mtime>=pdf.stat().st_mtime for x in pngs)
    names=['code/memory_charged_transfer.py','code/check_memory_polar_transfer.py',
           'code/verify_memory_transfer_closure.py','code/build_memory_charged_transfer.ps1',
           'code/record_memory_charged_transfer.py',
           'manuscript/susskind-memory-charged-transfer/main.tex',
           'output/pdf/de-sitter-susskind-memory-charged-transfer.pdf',
           'notes/susskind-memory-charged-transfer.md',
           'literature/susskind-memory-charged-transfer-sources.md',
           'journal/2026-11-07-susskind-memory-charged-transfer.md',
           'README.md','PLAN.md','NEXT.md']
    paths=[ROOT/x for x in names]+[OUT/x for x in data]+[OUT/'main.log']
    paths+=sorted((OUT/'diagnostics').glob('*'))
    result=dict(
        status='Local research progress; full research goal remains active',
        created_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        scientific_status='Actual adjacent-charge ordinary propagator approximated by same-original-disorder Gaussian OU smoothing with full-charge annealed L2 error O_S(1/N). Derivative/energy-weighted transfer, collective force defects, X0 covariance and bulk dictionary remain open.',
        git_commit_created=False,git_push_performed=False,archive_created=False,external_reply_sent=False,
        previous_evidence=parents,verified_scientific_inputs=verified,
        checks=dict(complete_CAR_systems=6,scientific_outputs=3,symbolic_identities=3,
                    exact_rational_source_variance_cases=230,max_matrix_errors=maxerrors),
        environment=dict(python=platform.python_version(),numpy=numpy.__version__,
                         scipy=scipy.__version__,sympy=sympy.__version__,mpmath=mpmath.__version__),
        parameters='Finite systems recorded in each JSON; C_upper=1 for CAR checks. Analytic DS result assumes lambda_N in a positive compact interval, fixed Js and compact real time; central finite estimates assume N>=16, p<=N/16, N/4<=m<=3N/4.',
        evidence_scopes={name:record['scope'] for name,record in data.items()},
        proof_scope='Analytic Sections 2-7 with independent read-only mathematical reviews. Small-system superoperator checks are not fits or substitutes for the uniform proof. The OU parameter alpha is not the energy-Mehler parameter q; no time differentiation of an unweighted error is asserted.',
        build=dict(command='powershell -NoProfile -ExecutionPolicy Bypass -File code/build_memory_charged_transfer.ps1',
                   pages=8,passes=2,compiler_header=log.splitlines()[0],
                   render_command='pdftoppm -r 95 -png output/pdf/de-sitter-susskind-memory-charged-transfer.pdf tmp/pdfs/susskind-memory-charged-transfer-build/page',
                   final_pages_visually_inspected=list(range(1,9)),
                   visual_result='Final eight pages have legible equations, clean margins and references, and no clipping or overlaps.'),
        diagnostics=['First compile failed because mathscr was not defined; source and log retained. Changed notation to mathcal; final two-pass compile clean.',
                     'One context-mismatched apply_patch was rejected without edits; re-read exact lines and applied corrected patch.',
                     'Root DOI open for Choi returned Internal Error; collaborator checked the author-uploaded original. Root directly read Kota and Chen-Liu formulas.',
                     'All three scientific drivers passed first execution; no old scientific JSON was overwritten.'],
        process_state='All scientific and build shell commands terminal; no live computation handles.',
        hashes={x.relative_to(ROOT).as_posix():sha(x) for x in paths},
        render_hashes={x.relative_to(ROOT).as_posix():sha(x) for x in pngs},
        build_pass_hashes={x.relative_to(ROOT).as_posix():sha(x) for x in passes})
    with dest.open('x',encoding='utf-8') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(dict(status='passed',pages=8,
                         previous_files_checked=sum(len(x['unchanged_scientific_files']) for x in parents.values()),
                         bound_stage_files=len(paths),max_errors=maxerrors)))
