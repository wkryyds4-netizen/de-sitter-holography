"""Bind the energy-regression and first memory-fluctuation stage locally."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib, json, re, subprocess, platform

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-fluctuations'
BUILD=ROOT/'tmp/pdfs/susskind-memory-fluctuations-build'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    target=OUT/'local-evidence.json'
    if target.exists():raise FileExistsError(target)
    parent=json.loads((ROOT/'results/susskind-memory-force/local-evidence.json').read_text())
    previous=[]
    for stage in [r['stage'] for r in parent['previous_stages']]+['susskind-memory-force']:
        path=ROOT/f'results/{stage}/local-evidence.json'
        old=json.loads(path.read_text());kept=[]
        for name,digest in old['hashes'].items():
            if name in ['README.md','PLAN.md','NEXT.md']:continue
            assert sha(ROOT/name)==digest,name
            kept.append(name)
        previous.append(dict(stage=stage,manifest_sha256=sha(path),unchanged_verified=kept))
    record=json.loads((OUT/'checks.json').read_text())
    assert record['status']=='passed'
    for path,digest in record['input_sha256'].items():assert sha(ROOT/path)==digest,path
    assert len(record['finite'])==3 and sum(r['exact_identities'] for r in record['finite'])==15
    assert len(record['symbolic']['conditional_moment_tests'])==28
    assert record['symbolic']['small_lambda_coefficients']==['4','32','768']
    assert len(record['curves'])==4
    assert sum('mehler_fine' in r for r in record['curves'])==2
    log=(OUT/'main.log').read_text(errors='replace')
    assert re.search(r'6\s+pages',log)
    assert not re.search(r'Overfull|undefined|Rerun to get|Fatal error',log)
    pdf=ROOT/'output/pdf/de-sitter-susskind-memory-fluctuations.pdf'
    assert sha(pdf)==sha(BUILD/'main.pdf')
    renders=[BUILD/f'page-{i}.png' for i in range(1,7)]
    assert all(p.exists() for p in renders)
    head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    assert head==parent['current_git_head']
    names=['code/memory_fluctuation_checks.py','code/build_memory_fluctuations.ps1',
           'code/record_memory_fluctuations.py','manuscript/susskind-memory-fluctuations/main.tex',
           'output/pdf/de-sitter-susskind-memory-fluctuations.pdf',
           'notes/susskind-memory-fluctuations.md','literature/susskind-memory-fluctuations-sources.md',
           'journal/2026-11-06-susskind-memory-fluctuations.md','README.md','PLAN.md','NEXT.md']
    files=[ROOT/n for n in names]+sorted(p for p in OUT.iterdir() if p.is_file())
    result=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),
        status='Local stage saved; no publication or global freeze',
        scientific_status='Energy-weighted mean, scalar-centering lower bound, exact first two conserved counterterms and complete spectral matrix of the memory coefficient pair. Finite-memory fluctuation convergence and bulk spatial dictionary remain open.',
        previous_stages=previous,current_git_head=head,git_commit_created=False,
        git_push_performed=False,archive_created=False,global_project_audit_performed=False,
        checks=dict(exact_finite_identities=15,conditional_polynomial_tests=28,
                    exact_laurent_coefficients=[4,32,768],lambda_curves=4,independent_mehler_cases=2),
        environment=dict(python=platform.python_version(),calculation_environment=record['environment']),
        build=dict(command='powershell -NoProfile -ExecutionPolicy Bypass -File code/build_memory_fluctuations.ps1',
                   compiler=log.splitlines()[0],passes=2,pages=6,
                   visual_review='All six final rendered pages inspected; no clipping, overlap or undefined citations',
                   render_command='pdftoppm -r 95 -png output/pdf/de-sitter-susskind-memory-fluctuations.pdf tmp/pdfs/susskind-memory-fluctuations-build/page'),
        retained_failures=['Initial Mehler tangent omitted K_(q^2) on the energy-descendant term. Source, failure explanation and numeric diagnostic retained. Factor restored without tolerance changes.'],
        error_scope=record['error_scope'],
        process_state='All computation, compilation and rendering processes terminal',
        hashes={p.relative_to(ROOT).as_posix():sha(p) for p in files},
        local_render_hashes={p.relative_to(ROOT).as_posix():sha(p) for p in renders},
        build_pass_hashes={f'pass{i}.txt':sha(BUILD/f'pass{i}.txt') for i in [1,2]})
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status='saved locally',files=len(files),pages=6,
                         prior_scientific_files=sum(len(r['unchanged_verified']) for r in previous),git_head=head)))


if __name__=='__main__':main()
