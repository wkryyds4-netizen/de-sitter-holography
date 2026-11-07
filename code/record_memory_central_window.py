"""Save checked local evidence without a Git mutation."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,platform,re,subprocess

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-central-window'
BUILD=ROOT/'tmp/pdfs/susskind-memory-central-window-build'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    target=OUT/'local-evidence.json'
    if target.exists():raise FileExistsError(target)
    records={}
    for name in ['central-moments','central-scale','central-filter-checks']:
        data=json.loads((OUT/(name+'.json')).read_text())
        assert data['status']=='passed'
        for n,h in data['input_sha256'].items():assert sha(ROOT/n)==h,n
        records[name]=data
    moments=records['central-moments'];scale=records['central-scale'];filters=records['central-filter-checks']
    assert len(moments['symbolic_rows'])==6 and moments['cancelled_negative_powers']==21
    assert len(moments['exact_full_distribution_checks'])==4 and len(moments['finite_cases'])==40
    assert len(scale['exact_variance_checks'])==30
    assert filters['exact_filter_identities']==60 and len(filters['records'])==30
    parentpath=ROOT/'results/susskind-memory-raw-spectrum/local-evidence.json'
    parent=json.loads(parentpath.read_text());unchanged=[]
    for n,h in parent['hashes'].items():
        if n in ['README.md','PLAN.md','NEXT.md']:continue
        assert sha(ROOT/n)==h,n
        unchanged.append(n)
    log=(OUT/'main.log').read_text(errors='replace')
    assert re.search(r'6\s+pages',log)
    assert not re.search(r'Overfull|undefined|Rerun to get|Fatal error',log)
    pdf=ROOT/'output/pdf/de-sitter-susskind-memory-central-window.pdf'
    assert sha(pdf)==sha(BUILD/'main.pdf')
    info=subprocess.check_output(['pdfinfo',str(pdf)],text=True,errors='replace')
    assert re.search(r'Pages:\s+6\b',info)
    renders=[BUILD/f'page-{i}.png' for i in range(1,7)]
    assert all(p.exists() for p in renders)
    head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    assert head==parent['current_git_head']
    names=['code/memory_central_window.py','code/memory_central_scale.py','code/check_memory_central_filter.py',
           'code/build_memory_central_window.ps1','code/record_memory_central_window.py',
           'manuscript/susskind-memory-central-window/main.tex',
           'output/pdf/de-sitter-susskind-memory-central-window.pdf',
           'notes/susskind-memory-central-window.md','literature/susskind-memory-central-window-sources.md',
           'journal/2026-11-06-susskind-memory-central-window.md','README.md','PLAN.md','NEXT.md']
    paths=[ROOT/n for n in names]+sorted(p for p in OUT.rglob('*') if p.is_file())
    result=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),status='Saved locally; no publication',
        scientific_status='Uniform compact-central-window sixth rank-product moment; explicit sector variance drift and finite conserved gH filter correction. Full memory Gram, charge-tail integration, finite-memory remainder and bulk dictionary remain open.',
        current_git_head=head,git_commit_created=False,git_push_performed=False,archive_created=False,
        parent_manifest_sha256=sha(parentpath),previous_stage_scientific_files_unchanged=unchanged,
        checks=dict(symbolic_moments=6,negative_power_cancellations=21,finite_full_distribution_moments=28,
                    central_moment_cases=40,variance_ratios=30,exact_filter_identities=60,filter_convergence_cases=30),
        environment=dict(python=platform.python_version(),calculations=[d['environment'] for d in records.values()]),
        build=dict(command='powershell -NoProfile -ExecutionPolicy Bypass -File code/build_memory_central_window.ps1',
                   compiler=log.splitlines()[0],passes_per_build=2,final_pages=6,
                   render_command='pdftoppm -r 95 -png output/pdf/de-sitter-susskind-memory-central-window.pdf tmp/pdfs/susskind-memory-central-window-build/page',
                   visual_review='All six final page PNGs viewed at 95 dpi; formula bounds, reference text and page boundaries clear.'),
        failures=[dict(stage='central symbolic recurrence',reason='Touchard dummy z shadowed polynomial-ring charge z; renamed tvar. Exact checks retained.',resolved=True),
                  dict(stage='scale symbolic equality',reason='Structural equality rejected algebraically identical factorization; changed to exact factor(difference)==0.',resolved=True),
                  dict(stage='first PDF build',reason='2.36888pt overflow in parameter list; saved source/log and moved list to display.',resolved=True)],
        evidence_scope='Rational analyticity and polynomial-ring identities prove compact-window moment estimates. Finite five-word checks verify exact filters. No uniform full memory Gram or finite-time remainder theorem is inferred.',
        process_state='All calculations, compilation and rendering terminal; no active calculation.',
        hashes={p.relative_to(ROOT).as_posix():sha(p) for p in paths},
        render_hashes={p.relative_to(ROOT).as_posix():sha(p) for p in renders},
        build_pass_hashes={f'pass{i}.txt':sha(BUILD/f'pass{i}.txt') for i in [1,2]})
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status='saved locally',pages=6,files=len(paths),
                         previous_scientific_files_unchanged=len(unchanged),git_head=head)))


if __name__=='__main__':main()
