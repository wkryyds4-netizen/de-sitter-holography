"""Record the completed cubic-cross stage locally, preserving previous science."""
from pathlib import Path
from datetime import datetime,timezone
import json,hashlib,re,subprocess,platform
import numpy,mpmath,sympy
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-cubic-cross'
BUILD=ROOT/'tmp/pdfs/susskind-cubic-cross-build'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    target=OUT/'local-evidence.json'
    if target.exists():raise FileExistsError(target)
    parent=json.loads((ROOT/'results/susskind-linear-kernel/local-evidence.json').read_text())
    stages=[x['stage'] for x in parent['previous_stages']]+['susskind-linear-kernel'];previous=[]
    for stage in stages:
        path=ROOT/f'results/{stage}/local-evidence.json';record=json.loads(path.read_text());kept=[]
        for name,h in record['hashes'].items():
            if name in ['README.md','PLAN.md','NEXT.md']:continue
            assert sha(ROOT/name)==h,name;kept.append(name)
        previous.append(dict(stage=stage,manifest_sha256=sha(path),unchanged_verified=kept))
    data={}
    for name in ['finite-cross','linear-channels','counterterm-limits','cross-orbits','reduction-checks']:
        record=json.loads((OUT/f'{name}.json').read_text());data[name]=record
        assert record['status'] in ['passed','computed']
        for path,h in record['input_sha256'].items():assert sha(ROOT/path)==h,path
    assert len(data['finite-cross']['records'])==3
    assert data['finite-cross']['records'][-1]['moments']['Z']['cross_M2']=='-2160'
    assert data['linear-channels']['checks']==32
    assert len(data['counterterm-limits']['records'])==15
    assert data['cross-orbits']['templates']==177 and data['cross-orbits']['words']==6234
    assert data['cross-orbits']['minor_checks']==37404
    assert len(data['reduction-checks']['records'])==3
    assert sum(r['exact_operator_identities'] for r in data['reduction-checks']['records'])==12
    failure=json.loads((OUT/'initial-reduction-check-failure.json').read_text())
    assert failure['status']=='failed' and (ROOT/failure['original_code']).exists()
    log=(OUT/'main.log').read_text(errors='replace')
    assert not re.search(r'Overfull|undefined|Rerun to get|Fatal error',log) and '(7 pages,' in log
    assert sha(ROOT/'output/pdf/de-sitter-susskind-cubic-cross.pdf')==sha(BUILD/'main.pdf')
    renders=[BUILD/f'final-page-{i}.png' for i in range(1,8)]
    assert all(p.exists() for p in renders)
    for i in range(1,7):assert sha(BUILD/f'page-{i}.png')==sha(BUILD/f'final-page-{i}.png')
    names=['code/cubic_cross_probe.py','code/cubic_cross_channels.py','code/run_cubic_cross_limits.py',
        'code/cubic_cross_orbits.py','code/check_cubic_cross.py','code/build_cubic_cross.ps1','code/record_cubic_cross.py',
        'manuscript/susskind-cubic-cross/main.tex','output/pdf/de-sitter-susskind-cubic-cross.pdf',
        'notes/susskind-cubic-cross.md','literature/susskind-cubic-cross-sources.md',
        'journal/2026-11-06-susskind-cubic-cross.md','README.md','PLAN.md','NEXT.md']
    files=[ROOT/n for n in names]+sorted(p for p in OUT.iterdir() if p.is_file())
    head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    assert head=='714d038bb605cfcffdc351e4dcc62aa95237d25f'
    record=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),
        status='Local stage saved; no publication or global freeze',
        scientific_status='Exact raw-cubic reduction, 1/p S/Z counterterms and F-Z leakage; all-order R-S/R-Z cancellation and vanishing R-Y two-energy cross measure. F-R and energy-resolved F-Y remain open.',
        previous_stages=previous,current_git_head=head,git_commit_created=False,git_push_performed=False,
        archive_created=False,global_project_audit_performed=False,
        environment=dict(python=platform.python_version(),numpy=numpy.__version__,mpmath=mpmath.__version__,sympy=sympy.__version__),
        checks=dict(finite_cross_systems=3,linear_channels=32,full_rank_cases=15,precisions=[80,120],
            orbit_templates=177,orbit_words=6234,local_minor_checks=37404,
            independent_rank_curvatures=3,exact_deterministic_matrix_identities=12,
            precision_scope='Full rank sums and exact polynomial/rational checks; no interval rounding certificate'),
        build=dict(command='powershell -NoProfile -ExecutionPolicy Bypass -File code/build_cubic_cross.ps1',
            compiler=log.splitlines()[0],passes=2,pages=7,
            visual_review='All seven pages inspected; final page 7 reinspected after adding validation evidence, pages 1-6 pixel-identical.'),
        retained_failure='Initial deterministic checker structural equality failed on unexpanded expressions; exact expansion fixed comparison, no tolerance or scientific change.',
        process_state='All calculation, compilation and rendering sessions terminal; resumed session52786 exit0',
        hashes={p.relative_to(ROOT).as_posix():sha(p) for p in files},
        local_render_hashes={p.relative_to(ROOT).as_posix():sha(p) for p in renders},
        build_pass_hashes={f'pass{i}.txt':sha(BUILD/f'pass{i}.txt') for i in [1,2]})
    target.write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(dict(status='saved locally',files=len(files),
        prior_scientific_files=sum(len(x['unchanged_verified']) for x in previous),pages=7,git_head=head)))
if __name__=='__main__':main()
