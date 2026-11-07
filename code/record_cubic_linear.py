"""Save local cubic-linear evidence and preserve prior scientific hashes."""
from pathlib import Path
from datetime import datetime,timezone
import json,hashlib,re,subprocess
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-cubic-linear'
BUILD=ROOT/'tmp/pdfs/susskind-cubic-linear-build'
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def main():
    target=OUT/'local-evidence.json'
    if target.exists():raise FileExistsError(target)
    stages=['susskind-singlet-moments','susskind-singlet-double-scale','susskind-singlet-kernel',
      'susskind-singlet-multimark','susskind-singlet-otoc','susskind-rindler-packet',
      'susskind-first-order','susskind-translations','susskind-antisymmetric',
      'susskind-two-matter','susskind-cubic-singlet','susskind-cubic-kernel']
    previous=[]
    for stage in stages:
        path=ROOT/f'results/{stage}/local-evidence.json';data=json.loads(path.read_text());kept=[]
        for name,h in data['hashes'].items():
            if name in ['README.md','PLAN.md','NEXT.md']:continue
            assert sha(ROOT/name)==h,name
            kept.append(name)
        previous.append(dict(stage=stage,manifest_sha256=sha(path),unchanged_verified=kept))
    for name in ['finite-probe.json','channels.json','validation.json','double-scale.json']:
        data=json.loads((OUT/name).read_text());assert data['status'] in ['passed','computed']
        for path,h in data['input_sha256'].items():assert sha(ROOT/path)==h,path
    check=json.loads((OUT/'validation.json').read_text())
    assert check['charged_covariance_checks']==31 and check['exact_weight_normalizations']==199
    assert check['witness']['normalized_commutator_squared']=='8933634/5'
    assert json.loads((OUT/'channels.json').read_text())['channel_checks']==42
    assert len(json.loads((OUT/'double-scale.json').read_text())['records'])==12
    log=(OUT/'main.log').read_text(errors='replace')
    assert not re.search(r'Overfull|undefined|Rerun to get|Fatal error',log) and '(6 pages,' in log
    assert sha(ROOT/'output/pdf/de-sitter-susskind-cubic-linear.pdf')==sha(BUILD/'main.pdf')
    renders=[BUILD/f'final-page-{i}.png' for i in range(1,7)]
    assert all(p.exists() for p in renders)
    for n in [1,2,4,5,6]:assert sha(BUILD/f'page-{n}.png')==sha(BUILD/f'final-page-{n}.png')
    files=[ROOT/p for p in [
      'code/cubic_linear_probe.py','code/cubic_linear_channels.py','code/cubic_linear_limit.py',
      'code/check_cubic_linear.py','code/build_cubic_linear.ps1','code/record_cubic_linear.py',
      'manuscript/susskind-cubic-linear/main.tex','output/pdf/de-sitter-susskind-cubic-linear.pdf',
      'notes/susskind-cubic-linear.md','literature/susskind-cubic-linear-sources.md',
      'journal/2026-11-06-susskind-cubic-linear.md','README.md','PLAN.md','NEXT.md']]
    files+=sorted(p for p in OUT.iterdir() if p.is_file())
    head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    assert head=='714d038bb605cfcffdc351e4dcc62aa95237d25f'
    data=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),
      status='Local stage saved; no publication or global freeze',
      scientific_status='Exact finite-N first-chaos coefficients, finite nonconservation witness and positive energy-orthogonal DS residual norm. Its dynamic spectral measure, mixed real-time correlators, bulk mass dictionary and full duality remain open.',
      previous_stages=previous,current_git_head=head,git_commit_created=False,git_push_performed=False,
      archive_created=False,global_project_audit_performed=False,
      checks=dict(channel_coefficients=42,charged_covariance=31,exact_weight_sets=199,symbolic_projection=3,
        double_scale_cases=12,precisions=[80,110],precision_comparison_scope='60-digit saved variance values'),
      build=dict(command='powershell -NoProfile -ExecutionPolicy Bypass -File code/build_cubic_linear.ps1',
        compiler=log.splitlines()[0],passes=2,pages=6,
        visual_review='All six pages inspected; final page 3 re-inspected after clarification, other pages pixel-identical.'),
      process_state='All calculation, compilation and rendering sessions terminal',
      hashes={p.relative_to(ROOT).as_posix():sha(p) for p in files},
      local_render_hashes={p.relative_to(ROOT).as_posix():sha(p) for p in renders},
      build_pass_hashes={f'pass{i}.txt':sha(BUILD/f'pass{i}.txt') for i in [1,2]})
    target.write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(dict(status='saved locally',files=len(files),prior_scientific_files=sum(len(p['unchanged_verified']) for p in previous),pages=6,git_head=head)))
if __name__=='__main__':main()
