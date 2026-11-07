"""Preserve the commutator-primary stage without Git publication."""
from pathlib import Path
from datetime import datetime,timezone
import json,hashlib,platform,re,subprocess
import numpy,mpmath,sympy,matplotlib
ROOT=Path(__file__).resolve().parent.parent;OUT=ROOT/'results/susskind-commutator-primary'
BUILD=ROOT/'tmp/pdfs/susskind-commutator-primary-build'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    target=OUT/'local-evidence.json'
    if target.exists():raise FileExistsError(target)
    parent=json.loads((ROOT/'results/susskind-mixed-cubic/local-evidence.json').read_text())
    stages=[r['stage'] for r in parent['previous_stages']]+['susskind-mixed-cubic'];previous=[]
    for stage in stages:
        path=ROOT/f'results/{stage}/local-evidence.json';r=json.loads(path.read_text());kept=[]
        for n,h in r['hashes'].items():
            if n in ['README.md','PLAN.md','NEXT.md']:continue
            assert sha(ROOT/n)==h,n;kept.append(n)
        previous.append(dict(stage=stage,manifest_sha256=sha(path),unchanged_verified=kept))
    data={}
    for n in ['exact','vertex','numerics']:
        r=json.loads((OUT/f'{n}.json').read_text());data[n]=r;assert r['status']=='passed'
        for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h,p
    exact=data['exact'];vertex=data['vertex'];num=data['numerics']
    assert len(exact['primary']['two_arc_cases'])==25 and len(exact['primary']['even_moments'])==3
    assert len(exact['covariance']['records'])==20 and exact['covariance']['laurent_cancellations']==60
    assert exact['primary']['fourth_difference_over_h2_limit']=='54' and vertex['three_arc_cases']==50
    assert len(num['records'])==3 and len(num['nonsymmetric_checks'])==4
    for r in num['records']:
        assert r['refinement_absolute_difference']<1e-10
        for p in r['points']:
            assert abs(p['coefficient']-p['commutator_squared']/r['lambda_value'])<1e-12
    log=(OUT/'main.log').read_text(errors='replace')
    assert not re.search(r'Overfull|undefined|Rerun to get|Fatal error',log) and '(6 pages,' in log
    assert '65.35826pt' in (OUT/'build-initial.log').read_text(errors='replace')
    pdf=ROOT/'output/pdf/de-sitter-susskind-commutator-primary.pdf';assert sha(pdf)==sha(BUILD/'main.pdf')
    renders=[BUILD/f'page-{i}.png' for i in range(1,7)];assert all(p.exists() for p in renders)
    head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip();assert head=='714d038bb605cfcffdc351e4dcc62aa95237d25f'
    names=['code/commutator_primary.py','code/check_commutator_vertex.py','code/run_commutator_primary.py',
      'code/build_commutator_primary.ps1','code/record_commutator_primary.py',
      'manuscript/susskind-commutator-primary/main.tex','output/pdf/de-sitter-susskind-commutator-primary.pdf',
      'notes/susskind-commutator-primary.md','literature/susskind-commutator-primary-sources.md',
      'journal/2026-11-06-susskind-commutator-primary.md','README.md','PLAN.md','NEXT.md']
    files=[ROOT/n for n in names]+sorted(p for p in OUT.iterdir() if p.is_file())
    record=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),status='Local stage saved; no publication or global freeze',
        scientific_status='Sourced fixed-time mixed growth. Exact k0 order-odd q^5 propagation and three-time vertex; composite fourth-moment excess. Full spatial/mass/interaction dictionary and de Sitter duality remain open.',
        previous_stages=previous,current_git_head=head,git_commit_created=False,git_push_performed=False,archive_created=False,
        global_project_audit_performed=False,
        checks=dict(two_arc_moments=25,three_arc_identities=50,even_moments=3,covariance_coefficients=20,
            laurent_cancellations=60,finite_lambda_groups=3,high_precision_transfer_points=4),
        environment=dict(python=platform.python_version(),numpy=numpy.__version__,mpmath=mpmath.__version__,sympy=sympy.__version__,matplotlib=matplotlib.__version__),
        build=dict(command='powershell -NoProfile -ExecutionPolicy Bypass -File code/build_commutator_primary.ps1',
          compiler=log.splitlines()[0],passes=2,pages=6,visual_review='All six final pages inspected',
          render_command='pdftoppm -r 95 -png output/pdf/de-sitter-susskind-commutator-primary.pdf tmp/pdfs/susskind-commutator-primary-build/page'),
        error_scope='Exact polynomial identities; analytic Jacobi truncation bounds. Floating-point controls empirical; sourced small-h expansion has no certified uniform numerical remainder here.',
        retained_failures=['Initial 65.35826pt typesetting overflow, fixed by splitting math list; log retained'],
        process_state='All computation, compilation and rendering processes terminal',
        hashes={p.relative_to(ROOT).as_posix():sha(p) for p in files},
        local_render_hashes={p.relative_to(ROOT).as_posix():sha(p) for p in renders},
        build_pass_hashes={f'pass{i}.txt':sha(BUILD/f'pass{i}.txt') for i in [1,2]})
    target.write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(dict(status='saved locally',files=len(files),prior_scientific_files=sum(len(p['unchanged_verified']) for p in previous),pages=6,git_head=head)))
if __name__=='__main__':main()
