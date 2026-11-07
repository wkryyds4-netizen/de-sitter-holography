"""Local central-Gram evidence, without committing or publishing."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,platform,re,subprocess
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-central-gram'
BUILD=ROOT/'tmp/pdfs/susskind-memory-central-gram-build'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    target=OUT/'local-evidence.json'
    if target.exists():raise FileExistsError(target)
    records={}
    for name in ['neutral-drift','polarization-car','gram-assembly']:
        data=json.loads((OUT/(name+'.json')).read_text());assert data['status']=='passed'
        for n,h in data['input_sha256'].items():assert sha(ROOT/n)==h,n
        records[name]=data
    assert records['neutral-drift']['exact_weight_ratios']==758
    assert len(records['neutral-drift']['records'])==36
    assert sum(r['ordered_gaussian_pairs'] for r in records['polarization-car']['rows'])==3217
    assert records['gram-assembly']['ordered_small_lambda_limit']=='1024*J_s**8'
    parentpath=ROOT/'results/susskind-memory-central-cross/local-evidence.json'
    parent=json.loads(parentpath.read_text());unchanged=[]
    for n,h in parent['hashes'].items():
        if n in ['README.md','PLAN.md','NEXT.md']:continue
        assert sha(ROOT/n)==h,n
        unchanged.append(n)
    log=(OUT/'main.log').read_text(errors='replace')
    assert re.search(r'6\s+pages',log)
    assert not re.search(r'Overfull|undefined|Rerun to get|Fatal error',log)
    pdf=ROOT/'output/pdf/de-sitter-susskind-memory-central-gram.pdf'
    assert sha(pdf)==sha(BUILD/'main.pdf')
    info=subprocess.check_output(['pdfinfo',str(pdf)],text=True,errors='replace')
    assert re.search(r'Pages:\s+6\b',info)
    renders=[BUILD/f'page-{i}.png' for i in range(1,7)]
    assert all(p.exists() for p in renders)
    head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    assert head==parent['current_git_head']
    names=['code/memory_central_polarization.py','code/check_memory_polarization_car.py',
           'code/memory_central_gram_assembly.py','code/build_memory_central_gram.ps1',
           'code/record_memory_central_gram.py','manuscript/susskind-memory-central-gram/main.tex',
           'output/pdf/de-sitter-susskind-memory-central-gram.pdf','notes/susskind-memory-central-gram.md',
           'literature/susskind-memory-central-gram-sources.md','journal/2026-11-06-susskind-memory-central-gram.md',
           'README.md','PLAN.md','NEXT.md']
    paths=[ROOT/n for n in names]+sorted(p for p in OUT.rglob('*') if p.is_file())
    result=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),status='Saved locally; no publication',
        scientific_status='Even N compact central window: neutral drift O(1/(N p)), polarization proof of S--O, full four-field equal-time Gram, and projected W2 core variance. No assumption of the optional exact charged Casimir candidate. Energy-inserted spectra, full charge integration, finite-memory remainders and bulk dictionary remain open.',
        current_git_head=head,git_commit_created=False,git_push_performed=False,archive_created=False,
        parent_manifest_sha256=sha(parentpath),previous_stage_scientific_files_unchanged=unchanged,
        checks=dict(exact_weight_ratios=758,finite_rational_moments=4,all_rank_cases=36,
                    precisions=[80,120],CAR_systems=3,ordered_CAR_pairs=3217,rank_blocks=18),
        environment=dict(python=platform.python_version(),calculations=[d['environment'] for d in records.values()]),
        build=dict(command='powershell -NoProfile -ExecutionPolicy Bypass -File code/build_memory_central_gram.ps1',
                   compiler=log.splitlines()[0],passes_in_final_build=2,final_pages=6,
                   render_command='pdftoppm -r 95 -png output/pdf/de-sitter-susskind-memory-central-gram.pdf tmp/pdfs/susskind-memory-central-gram-build/page',
                   visual_review='All six final page PNGs viewed at 95 dpi; formulas, citations and page boundaries checked.'),
        failures=[dict(stage='symbolic score assertion',reason='Structural equality rejected algebraically identical factorizations; exact factored difference is zero. Preserved source/diagnostic, replaced structural comparison.',resolved=True),
                  dict(stage='first LaTeX build',reason='One unescaped fraction command became a form-feed character. Preserved source/log; escaped command correctly and completed two passes.',resolved=True)],
        evidence_scope='Analytic uniform neutral bounds and exact polarization establish the limiting cross; finite independent CAR and arbitrary-precision full sums check the formulas. The earlier optional finite charged identity remains unproved and is not used here.',
        process_state='All calculations, compilation and rendering terminal. Session 89893 confirmed exit 0.',
        hashes={p.relative_to(ROOT).as_posix():sha(p) for p in paths},
        render_hashes={p.relative_to(ROOT).as_posix():sha(p) for p in renders},
        build_pass_hashes={f'pass{i}.txt':sha(BUILD/f'pass{i}.txt') for i in [1,2]})
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status='saved locally',pages=6,files=len(paths),
                         previous_scientific_files_unchanged=len(unchanged),git_head=head)))

if __name__=='__main__':main()
