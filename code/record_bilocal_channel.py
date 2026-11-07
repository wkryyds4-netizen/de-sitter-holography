"""Bind the bilocal-channel stage and verify retained scientific history."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib,json,re,subprocess,platform
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-bilocal-channel'
BUILD=ROOT/'tmp/pdfs/susskind-bilocal-channel-build'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    target=OUT/'local-evidence.json'
    if target.exists():raise FileExistsError(target)
    parent=json.loads((ROOT/'results/susskind-bilocal-interface/local-evidence.json').read_text())
    previous=[]
    for stage in [r['stage'] for r in parent['previous_stages']]+['susskind-bilocal-interface']:
        path=ROOT/f'results/{stage}/local-evidence.json';r=json.loads(path.read_text());kept=[]
        for n,h in r['hashes'].items():
            if n in ['README.md','PLAN.md','NEXT.md']:continue
            assert sha(ROOT/n)==h,n
            kept.append(n)
        previous.append(dict(stage=stage,manifest_sha256=sha(path),unchanged_verified=kept))
    records={}
    for n in ['checks','projection-checks']:
        r=json.loads((OUT/f'{n}.json').read_text());assert r['status']=='passed'
        for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h,p
        records[n]=r
    c=records['checks'];p=records['projection-checks']
    assert len(c['deterministic'])==4 and len(c['covariance'])==5
    assert sum(x['exact_identities'] for x in c['deterministic'])==20
    assert len(c['double_scaling'])==18 and len(p['records'])==2
    log=(OUT/'main.log').read_text(errors='replace')
    assert re.search(r'6\s+pages',log) and not re.search(r'Overfull|undefined|Rerun to get|Fatal error',log)
    pdf=ROOT/'output/pdf/de-sitter-susskind-bilocal-channel.pdf'
    assert sha(pdf)==sha(BUILD/'main.pdf')
    renders=[BUILD/f'page-{i}.png' for i in range(1,7)]
    assert all(x.exists() for x in renders)
    head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    assert head==parent['current_git_head']
    names=['code/bilocal_channel.py','code/check_bilocal_projection.py',
      'code/build_bilocal_channel.ps1','code/record_bilocal_channel.py',
      'manuscript/susskind-bilocal-channel/main.tex',
      'output/pdf/de-sitter-susskind-bilocal-channel.pdf',
      'notes/susskind-bilocal-channel.md','literature/susskind-bilocal-channel-sources.md',
      'journal/2026-11-06-susskind-bilocal-channel.md','README.md','PLAN.md','NEXT.md']
    files=[ROOT/n for n in names]+sorted(x for x in OUT.iterdir() if x.is_file())
    result=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),
      status='Local stage saved; no publication or global freeze',
      scientific_status='Exact bilocal dilation and memory equation, initial singlet fluctuation identification, uniform compact-time trivial limit of ordinarily normalized bilocal. Full rescaled fluctuation kernel and holographic spatial map remain open.',
      previous_stages=previous,current_git_head=head,git_commit_created=False,
      git_push_performed=False,archive_created=False,global_project_audit_performed=False,
      checks=dict(integer_matrix_identities=20,exact_covariance_cases=5,
        high_precision_rows=18,finite_time_bounds=12,independent_projected_systems=2),
      environment=dict(python=platform.python_version(),calculation_environments=[r['environment'] for r in records.values()]),
      build=dict(command='powershell -NoProfile -ExecutionPolicy Bypass -File code/build_bilocal_channel.ps1',
        compiler=log.splitlines()[0],passes=2,pages=6,
        visual_review='All six final rendered pages inspected; no clipping, overlaps or undefined citations',
        render_command='pdftoppm -r 95 -png output/pdf/de-sitter-susskind-bilocal-channel.pdf tmp/pdfs/susskind-bilocal-channel-build/page'),
      retained_failures=['Initial LaTeX omitted the tr macro; original source and failure log preserved. Macro added; final build passed.'],
      error_scope={n:r['error_scope'] for n,r in records.items()},
      process_state='All computation, compilation and rendering processes terminal',
      hashes={x.relative_to(ROOT).as_posix():sha(x) for x in files},
      local_render_hashes={x.relative_to(ROOT).as_posix():sha(x) for x in renders},
      build_pass_hashes={f'pass{i}.txt':sha(BUILD/f'pass{i}.txt') for i in [1,2]})
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status='saved locally',files=len(files),pages=6,
      prior_scientific_files=sum(len(x['unchanged_verified']) for x in previous),git_head=head)))

if __name__=='__main__':main()
