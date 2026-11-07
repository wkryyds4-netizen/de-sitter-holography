"""Bind locally saved bilocal-interface computations, PDF and predecessor hashes."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib, json, re, subprocess, platform
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-bilocal-interface'
BUILD=ROOT/'tmp/pdfs/susskind-bilocal-interface-build'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    target=OUT/'local-evidence.json'
    if target.exists():raise FileExistsError(target)
    parent=json.loads((ROOT/'results/susskind-bulk-target/local-evidence.json').read_text())
    previous=[]
    for stage in [r['stage'] for r in parent['previous_stages']]+['susskind-bulk-target']:
        path=ROOT/f'results/{stage}/local-evidence.json'
        r=json.loads(path.read_text());kept=[]
        for n,h in r['hashes'].items():
            if n in ['README.md','PLAN.md','NEXT.md']:continue
            assert sha(ROOT/n)==h,n
            kept.append(n)
        previous.append(dict(stage=stage,manifest_sha256=sha(path),unchanged_verified=kept))
    records={}
    for n in ['micro','projection']:
        r=json.loads((OUT/f'{n}.json').read_text());assert r['status']=='passed'
        for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h,p
        records[n]=r
    m=records['micro'];p=records['projection']
    assert len(m['algebra'])==5 and len(m['combinatorial'])==77
    assert sum(r['monomials'] for r in m['direct_ensemble'])==54
    assert len(m['double_scaling'])==15 and len(p['sources'])==4
    assert all(r['quadrature_max_absolute_difference']<2e-12 for r in p['sources'])
    log=(OUT/'main.log').read_text(errors='replace')
    assert '(6 pages,' in log and not re.search(r'Overfull|undefined|Rerun to get|Fatal error',log)
    pdf=ROOT/'output/pdf/de-sitter-susskind-bilocal-interface.pdf'
    assert sha(pdf)==sha(BUILD/'main.pdf')
    renders=[BUILD/f'page-{i}.png' for i in range(1,7)]
    assert all(x.exists() for x in renders)
    head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    assert head==parent['current_git_head']
    names=['code/temporal_bilocal_interface.py','code/bilocal_meson_projection.py',
       'code/build_bilocal_interface.ps1','code/record_bilocal_interface.py',
       'manuscript/susskind-bilocal-interface/main.tex',
       'output/pdf/de-sitter-susskind-bilocal-interface.pdf',
       'notes/susskind-bilocal-interface.md','literature/susskind-bilocal-interface-sources.md',
       'journal/2026-11-06-susskind-bilocal-interface.md','README.md','PLAN.md','NEXT.md']
    files=[ROOT/n for n in names]+sorted(x for x in OUT.iterdir() if x.is_file())
    result=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),
       status='Local stage saved; no publication or global freeze',
       scientific_status='Known bulk projection applied to explicit sources; exact temporal-CAR degree weights and DS zero-time curvature. Full singlet projection and microscopic spatial map remain open.',
       previous_stages=previous,current_git_head=head,git_commit_created=False,
       git_push_performed=False,archive_created=False,global_project_audit_performed=False,
       checks=dict(exact_hamiltonians=5,combinatorial_cases=77,ensemble_monomials=54,
          high_precision_rows=15,bulk_sources=4),
       environment=dict(python=platform.python_version(),calculation_environments=[r['environment'] for r in records.values()]),
       build=dict(command='powershell -NoProfile -ExecutionPolicy Bypass -File code/build_bilocal_interface.ps1',
          compiler=log.splitlines()[0],passes=2,pages=6,
          visual_review='All six final rendered pages inspected; no clipping or collisions',
          render_command='pdftoppm -r 95 -png output/pdf/de-sitter-susskind-bilocal-interface.pdf tmp/pdfs/susskind-bilocal-interface-build/page',
          nonfatal_warning='One Underfull hbox in the bibliography; visually checked'),
       retained_failures=[],error_scope={n:r['error_scope'] for n,r in records.items()},
       process_state='All computation, compilation and rendering processes terminal',
       hashes={x.relative_to(ROOT).as_posix():sha(x) for x in files},
       local_render_hashes={x.relative_to(ROOT).as_posix():sha(x) for x in renders},
       build_pass_hashes={f'pass{i}.txt':sha(BUILD/f'pass{i}.txt') for i in [1,2]})
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status='saved locally',files=len(files),pages=6,
       prior_scientific_files=sum(len(x['unchanged_verified']) for x in previous),git_head=head)))

if __name__=='__main__':main()
