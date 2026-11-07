"""Save local cubic-kernel evidence and verify preceding scientific artifacts."""
from pathlib import Path
from datetime import datetime,timezone
import json,hashlib,re,subprocess
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-cubic-kernel'
BUILD=ROOT/'tmp/pdfs/susskind-cubic-kernel-build'
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def main():
    target=OUT/'local-evidence.json'
    if target.exists():raise FileExistsError(target)
    stages=['susskind-singlet-moments','susskind-singlet-double-scale','susskind-singlet-kernel',
      'susskind-singlet-multimark','susskind-singlet-otoc','susskind-rindler-packet',
      'susskind-first-order','susskind-translations','susskind-antisymmetric',
      'susskind-two-matter','susskind-cubic-singlet']
    previous=[]
    for stage in stages:
        path=ROOT/f'results/{stage}/local-evidence.json';data=json.loads(path.read_text());kept=[]
        for name,h in data['hashes'].items():
            if name in ['README.md','PLAN.md','NEXT.md']:continue
            assert sha(ROOT/name)==h,name
            kept.append(name)
        previous.append(dict(stage=stage,manifest_sha256=sha(path),unchanged_verified=kept))
    for name in ['orbit-certificate.json','independent-orbits.json','kernel.json','fourier-check.json']:
        data=json.loads((OUT/name).read_text());assert data['status']=='passed'
        for path,h in data['input_sha256'].items():assert sha(ROOT/path)==h,path
    orbit=json.loads((OUT/'orbit-certificate.json').read_text())
    check=json.loads((OUT/'independent-orbits.json').read_text())
    assert orbit['templates']==check['templates']==177
    assert orbit['orbit_words']==check['matching_words']==6234 and not orbit['failures']
    assert check['cancelled_templates']==176 and check['minor_calls']==56106
    kernel=json.loads((OUT/'kernel.json').read_text())
    assert len(kernel['records'])==5 and len(kernel['independent_mehler'])==2
    assert max(r['analytic_truncation_bound'] for r in kernel['records'])<1e-13
    fourier=json.loads((OUT/'fourier-check.json').read_text())
    assert fourier['checks_per_precision']==7
    log=(OUT/'main.log').read_text(errors='replace')
    assert not re.search(r'Overfull|undefined|Rerun to get|Fatal error',log)
    assert '(6 pages,' in log
    assert sha(ROOT/'output/pdf/de-sitter-susskind-cubic-kernel.pdf')==sha(BUILD/'main.pdf')
    renders=[BUILD/f'page-{i}.png' for i in range(1,7)]
    assert all(p.exists() for p in renders)
    files=[ROOT/p for p in [
      'code/cubic_kernel_orbits.py','code/check_cubic_kernel_orbits.py','code/run_cubic_kernel.py',
      'code/check_cubic_fourier.py','code/build_cubic_kernel.ps1','code/record_cubic_kernel.py',
      'manuscript/susskind-cubic-kernel/main.tex','output/pdf/de-sitter-susskind-cubic-kernel.pdf',
      'notes/susskind-cubic-kernel.md','literature/susskind-cubic-kernel-sources.md',
      'journal/2026-11-06-susskind-cubic-kernel.md','README.md','PLAN.md','NEXT.md']]
    files+=sorted(p for p in OUT.iterdir() if p.is_file())
    head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    assert head=='714d038bb605cfcffdc351e4dcc62aa95237d25f'
    data=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),
      status='Local stage saved; no publication or global freeze',
      scientific_status='All-order two-insertion q^3 reduction with exhaustive finite polynomial certificate; full fixed-time kernel and small-lambda sech^6 limit. Raw linear chaos, mixed higher correlators, physical mass dictionary and full duality remain open.',
      previous_stages=previous,current_git_head=head,git_commit_created=False,git_push_performed=False,
      archive_created=False,global_project_audit_performed=False,
      checks=dict(endpoint_templates=177,matching_words=6234,independent_minor_calls=56106,
        original_minor_cross_checks=3186,kernel_parameter_sets=5,mehler_parameter_sets=2,
        fourier_checks_per_precision=7,fourier_precisions=[50,75]),
      build=dict(command='powershell -NoProfile -ExecutionPolicy Bypass -File code/build_cubic_kernel.ps1',
        compiler=log.splitlines()[0],passes=2,pages=6,
        render_command='pdftoppm -scale-to 1100 -png output/pdf/de-sitter-susskind-cubic-kernel.pdf tmp/pdfs/susskind-cubic-kernel-build/page',
        visual_review='All six final pages visually inspected; equations, plot, table, bibliography and page edges clear.'),
      process_state='All calculation, compilation and rendering sessions terminal',
      hashes={p.relative_to(ROOT).as_posix():sha(p) for p in files},
      local_render_hashes={p.relative_to(ROOT).as_posix():sha(p) for p in renders},
      build_pass_hashes={f'pass{i}.txt':sha(BUILD/f'pass{i}.txt') for i in [1,2]})
    target.write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(dict(status='saved locally',files=len(files),prior_scientific_files=sum(len(p['unchanged_verified']) for p in previous),pages=6,git_head=head)))
if __name__=='__main__':main()
