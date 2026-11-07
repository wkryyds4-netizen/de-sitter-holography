"""Preserve the quartic-memory stage, all diagnostics and prior scientific files."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib,json,re,subprocess,platform

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-quartic'
BUILD=ROOT/'tmp/pdfs/susskind-memory-quartic-build'


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    target=OUT/'local-evidence.json'
    if target.exists():raise FileExistsError(target)
    parent=json.loads((ROOT/'results/susskind-memory-fluctuations/local-evidence.json').read_text())
    previous=[]
    for stage in [r['stage'] for r in parent['previous_stages']]+['susskind-memory-fluctuations']:
        path=ROOT/f'results/{stage}/local-evidence.json';old=json.loads(path.read_text());kept=[]
        for name,digest in old['hashes'].items():
            if name in ['README.md','PLAN.md','NEXT.md']:continue
            assert sha(ROOT/name)==digest,name
            kept.append(name)
        previous.append(dict(stage=stage,manifest_sha256=sha(path),unchanged_verified=kept))
    records={}
    for name in ['quartic','orbit-certificate','independent-car','cross-kernel']:
        record=json.loads((OUT/f'{name}.json').read_text());assert record['status']=='passed'
        for path,digest in record['input_sha256'].items():assert sha(ROOT/path)==digest,path
        records[name]=record
    o=records['orbit-certificate'];c=records['cross-kernel'];f=records['independent-car']
    assert o['templates']==807 and o['orbit_words']==452760 and not o['counterexamples']
    assert [r['templates'] for r in c['cross_certificates']]==[128,354]
    assert not any(r['counterexamples'] for r in c['cross_certificates'])
    assert c['surviving_orbit_matches_gram'] and c['small_lambda_variance_coefficient']=='768/5'
    assert len(f['finite'])==3 and sum(r['identities'] for r in f['finite'])==12
    assert f['local_words']['comparisons']==1600 and f['local_words']['nonzero']==239 and f['local_words']['negative']==8
    assert f['finite'][1]['zero_memory_control'] and f['finite'][1]['norm_squared'][2]=='0'
    log=(OUT/'main.log').read_text(errors='replace')
    assert re.search(r'6\s+pages',log)
    assert not re.search(r'Overfull|undefined|Rerun to get|Fatal error',log)
    pdf=ROOT/'output/pdf/de-sitter-susskind-memory-quartic.pdf'
    assert sha(pdf)==sha(BUILD/'main.pdf')
    renders=[BUILD/f'page-{i}.png' for i in range(1,7)]
    assert all(p.exists() for p in renders)
    head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    assert head==parent['current_git_head']
    names=['code/memory_quartic.py','code/quartic_orbits.py','code/check_memory_quartic_car.py',
           'code/quartic_cross_kernel.py','code/build_memory_quartic.ps1','code/record_memory_quartic.py',
           'manuscript/susskind-memory-quartic/main.tex','output/pdf/de-sitter-susskind-memory-quartic.pdf',
           'notes/susskind-memory-quartic.md','literature/susskind-memory-quartic-sources.md',
           'journal/2026-11-06-susskind-memory-quartic.md','README.md','PLAN.md','NEXT.md']
    files=[ROOT/n for n in names]+sorted(p for p in OUT.iterdir() if p.is_file())
    result=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),
        status='Local stage saved; no publication or global freeze',
        scientific_status='Exact second memory coefficient and return-term separation; quartic-chaos residual complete q^4 two-point kernel and zero energy-resolved Y/R crosses. Raw quadratic chaos, finite-memory covariance and bulk dictionary remain open.',
        previous_stages=previous,current_git_head=head,git_commit_created=False,
        git_push_performed=False,archive_created=False,global_project_audit_performed=False,
        checks=dict(endpoint_templates=807,endpoint_words=452760,cross_templates=482,
                    independent_finite_identities=12,finite_car_comparisons=1600,nonzero_car_values=239,negative_car_values=8,
                    kernel_curves=4,independent_quadrature_cases=2),
        environment=dict(python=platform.python_version(),calculation_environments=[r['environment'] for r in records.values()]),
        build=dict(command='powershell -NoProfile -ExecutionPolicy Bypass -File code/build_memory_quartic.ps1',
                   compiler=log.splitlines()[0],passes=2,pages=6,
                   visual_review='All six final pages inspected at 95 dpi; no clipping, overlap or undefined citations',
                   render_command='pdftoppm -r 95 -png output/pdf/de-sitter-susskind-memory-quartic.pdf tmp/pdfs/susskind-memory-quartic-build/page'),
        retained_issues=['Initial rounded Gaussian examples included two zero p=2 Hamiltonians; retained and superseded for coverage by independently verified nonzero integer couplings.',
                         'Unexpanded symbolic positivity check failed; subsequent slow nested-expression run deliberately interrupted and source preserved; per-step exact expansion then passed.',
                         'Parser rejected reserved token lambda; renamed parser token without changing saved expression.',
                         'LaTeX missing mathscr macro and later 3.16pt overfull paragraph; original source/logs retained, final six-page build passed.'],
        error_scope='All-order endpoint argument supported by complete exact local templates, not inferred from low moments. Jacobi truncation analytic; floating-point and quadrature refinement not interval-certified. Pure fourth chaos is not the complete raw second memory coefficient.',
        process_state='All computation, compilation and rendering terminal; deliberately terminated session78527 recorded; replacement session14677 exited zero',
        hashes={p.relative_to(ROOT).as_posix():sha(p) for p in files},
        local_render_hashes={p.relative_to(ROOT).as_posix():sha(p) for p in renders},
        build_pass_hashes={f'pass{i}.txt':sha(BUILD/f'pass{i}.txt') for i in [1,2]})
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status='saved locally',files=len(files),pages=6,
                         prior_scientific_files=sum(len(r['unchanged_verified']) for r in previous),git_head=head)))


if __name__=='__main__':main()
