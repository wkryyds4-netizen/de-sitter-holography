"""Bind the bulk-target stage, retained failures and previous scientific records."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib, json, re, subprocess, platform
import numpy as np
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-bulk-target'
BUILD=ROOT/'tmp/pdfs/susskind-bulk-target-build'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    target=OUT/'local-evidence.json'
    if target.exists():raise FileExistsError(target)
    pp=ROOT/'results/susskind-ordering-metric/local-evidence.json'
    parent=json.loads(pp.read_text());previous=[]
    for stage in [r['stage'] for r in parent['previous_stages']]+['susskind-ordering-metric']:
        path=ROOT/f'results/{stage}/local-evidence.json';r=json.loads(path.read_text());kept=[]
        for n,h in r['hashes'].items():
            if n in ['README.md','PLAN.md','NEXT.md']:continue
            assert sha(ROOT/n)==h,n;kept.append(n)
        previous.append(dict(stage=stage,manifest_sha256=sha(path),unchanged_verified=kept))
    records={}
    for n in ['spectrum','checks','mass-bound']:
        r=json.loads((OUT/f'{n}.json').read_text());assert r['status']=='passed'
        for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h,p
        records[n]=r
    spec=records['spectrum'];checks=records['checks']
    assert len(spec['chebyshev'])==4 and len(spec['fredholm'])==6
    assert len(checks['finite_part'])==12 and len(checks['angular_gram'])==36
    assert len(checks['inverse_moments'])==4 and len(records['mass-bound']['records'])==8
    assert float(checks['cross_representation_max_relative_error'])<7e-13
    with np.load(OUT/'wavefunctions.npz') as w:
        assert w['phi'].shape==(999,6) and w['coefficients'].shape==(256,6)
        assert np.allclose(w['twice_lambda'],spec['chebyshev'][-1]['twice_lambda'],rtol=0,atol=0)
        assert np.isfinite(w['phi']).all()
    log=(OUT/'main.log').read_text(errors='replace')
    assert '(5 pages,' in log and not re.search(r'Overfull|undefined|Rerun to get|Fatal error',log)
    pdf=ROOT/'output/pdf/de-sitter-susskind-bulk-target.pdf';assert sha(pdf)==sha(BUILD/'main.pdf')
    renders=[BUILD/f'page-{i}.png' for i in range(1,6)];assert all(p.exists() for p in renders)
    head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    assert head==parent['current_git_head']
    names=['code/thooft_bulk_target.py','code/check_thooft_bulk_target.py','code/check_bulk_mass_bound.py',
      'code/build_bulk_target.ps1','code/record_bulk_target.py','manuscript/susskind-bulk-target/main.tex',
      'output/pdf/de-sitter-susskind-bulk-target.pdf','notes/susskind-bulk-target.md',
      'literature/susskind-bulk-target-sources.md','journal/2026-11-06-susskind-bulk-target.md',
      'README.md','PLAN.md','NEXT.md']
    files=[ROOT/n for n in names]+sorted(p for p in OUT.iterdir() if p.is_file())
    result=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),status='Local stage saved; no publication or global freeze',
        scientific_status='Reproduction of the existing planar bulk meson target, analytic finite-parameter form bound, and conditional spatial-response mass test. Microscopic reconstruction and full duality remain open.',
        previous_stages=previous,current_git_head=head,git_commit_created=False,git_push_performed=False,
        archive_created=False,global_project_audit_performed=False,
        checks=dict(chebyshev_sizes=4,fredholm_domains_and_grids=6,levels=6,finite_part=12,gram=36,inverse_moments=4,finite_basis_bounds=8),
        environment=dict(python=platform.python_version(),calculation_environments=[r['environment'] for r in records.values()]),
        build=dict(command='powershell -NoProfile -ExecutionPolicy Bypass -File code/build_bulk_target.ps1',
          compiler=log.splitlines()[0],passes=2,pages=5,visual_review='All five final rendered pages inspected',
          render_command='pdftoppm -r 95 -png output/pdf/de-sitter-susskind-bulk-target.pdf tmp/pdfs/susskind-bulk-target-build/page'),
        retained_failures=['L32 Fredholm tail exceeded 2e-9 threshold; expanded to L40',
                           'mpmath chebyu hypergeometric evaluation failed; finite polynomial recurrence used'],
        error_scope=spec['error_scope'],process_state='All computation, compilation and rendering processes terminal',
        hashes={p.relative_to(ROOT).as_posix():sha(p) for p in files},
        local_render_hashes={p.relative_to(ROOT).as_posix():sha(p) for p in renders},
        build_pass_hashes={f'pass{i}.txt':sha(BUILD/f'pass{i}.txt') for i in [1,2]})
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status='saved locally',files=len(files),pages=5,
        prior_scientific_files=sum(len(p['unchanged_verified']) for p in previous),git_head=head)))

if __name__=='__main__':main()
