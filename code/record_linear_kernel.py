"""Save the covariance-bridge/full-kernel stage without a Git operation."""
from pathlib import Path
from datetime import datetime,timezone
import json,hashlib,re,subprocess
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-linear-kernel'
BUILD=ROOT/'tmp/pdfs/susskind-linear-kernel-build'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    target=OUT/'local-evidence.json'
    if target.exists():raise FileExistsError(target)
    parent=json.loads((ROOT/'results/susskind-cubic-dynamics/local-evidence.json').read_text())
    stages=[x['stage'] for x in parent['previous_stages']]+['susskind-cubic-dynamics'];previous=[]
    for stage in stages:
        path=ROOT/f'results/{stage}/local-evidence.json';record=json.loads(path.read_text());kept=[]
        for name,h in record['hashes'].items():
            if name in ['README.md','PLAN.md','NEXT.md']:continue
            assert sha(ROOT/name)==h,name;kept.append(name)
        previous.append(dict(stage=stage,manifest_sha256=sha(path),unchanged_verified=kept))
    for name in ['covariance-bridge.json','sixth-wick.json','joint-quartic.json','kernel.json']:
        record=json.loads((OUT/name).read_text());assert record['status']=='passed'
        for path,h in record['input_sha256'].items():assert sha(ROOT/path)==h,path
    bridge=json.loads((OUT/'covariance-bridge.json').read_text())
    assert bridge['exact_ratio_checks']==2286 and len(bridge['records'])==20
    sixth=json.loads((OUT/'sixth-wick.json').read_text())
    assert len(sixth['finite_records'])==3 and sixth['finite_records'][-1]['M4']=='314280'
    assert len(json.loads((OUT/'joint-quartic.json').read_text())['records'])==9
    kernel=json.loads((OUT/'kernel.json').read_text())
    assert len(kernel['records'])==5 and len(kernel['independent_mehler'])==2 and len(kernel['fourier'])==4
    failure=json.loads((OUT/'initial-plot-failure.json').read_text())
    assert failure['status']=='failed' and failure['scientific_results_written']==False
    log=(OUT/'main.log').read_text(errors='replace')
    assert not re.search(r'Overfull|undefined|Rerun to get|Fatal error',log) and '(7 pages,' in log
    assert sha(ROOT/'output/pdf/de-sitter-susskind-linear-kernel.pdf')==sha(BUILD/'main.pdf')
    renders=[BUILD/f'final-page-{i}.png' for i in range(1,8)]
    assert all(p.exists() for p in renders)
    for i in range(1,7):assert sha(BUILD/f'page-{i}.png')==sha(BUILD/f'final-page-{i}.png')
    names=['code/cubic_covariance_bridge.py','code/cubic_sixth_wick.py','code/check_linear_joint.py',
      'code/run_linear_kernel.py','code/build_linear_kernel.ps1','code/record_linear_kernel.py',
      'manuscript/susskind-linear-kernel/main.tex','output/pdf/de-sitter-susskind-linear-kernel.pdf',
      'notes/susskind-linear-kernel.md','literature/susskind-linear-kernel-sources.md',
      'journal/2026-11-06-susskind-linear-kernel.md','README.md','PLAN.md','NEXT.md']
    files=[ROOT/n for n in names]+sorted(p for p in OUT.iterdir() if p.is_file())
    head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    assert head=='714d038bb605cfcffdc351e4dcc62aa95237d25f'
    record=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),status='Local stage saved; no publication or global freeze',
      scientific_status='Rank-weight covariance bridge proves the all-fixed-degree orthogonal q-Gaussian (H,F) law using the existing independent embedded-ensemble chord rule. Full r=q Mehler kernel, no DS static atom and sech^2 small-lambda limit. Mixed F/R/S/Z spectra, bulk dictionary and full duality remain open.',
      previous_stages=previous,current_git_head=head,git_commit_created=False,git_push_performed=False,
      archive_created=False,global_project_audit_performed=False,
      checks=dict(exact_ratio=2286,full_rank_covariance_cases=20,finite_sixth_order_sectors=3,
        full_rank_joint_quartics=9,kernel_parameters=5,independent_mehler_cases=2,fourier_cases=4,
        precision_scope='80/120-digit rank sums, 65-digit stored comparisons; 50/80-digit Fourier; no interval rounding certification'),
      build=dict(command='powershell -NoProfile -ExecutionPolicy Bypass -File code/build_linear_kernel.ps1',
        compiler=log.splitlines()[0],passes=2,pages=7,
        visual_review='All seven final pages inspected. Only page 7 changed after adding joint-quartic evidence; pages 1-6 pixel-identical to first render.'),
      retained_failure='initial-plot-failure.json: corrected axis-label mathtext, unchanged scientific tolerances',
      process_state='All calculation, compilation and rendering sessions terminal',
      hashes={p.relative_to(ROOT).as_posix():sha(p) for p in files},
      local_render_hashes={p.relative_to(ROOT).as_posix():sha(p) for p in renders},
      build_pass_hashes={f'pass{i}.txt':sha(BUILD/f'pass{i}.txt') for i in [1,2]})
    target.write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(dict(status='saved locally',files=len(files),prior_scientific_files=sum(len(x['unchanged_verified']) for x in previous),pages=7,git_head=head)))
if __name__=='__main__':main()
