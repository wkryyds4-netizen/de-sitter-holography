"""Local central coherent-norm evidence; no Git writes."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,platform,re,subprocess
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-central-norm'
BUILD=ROOT/'tmp/pdfs/susskind-memory-central-norm-build'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    target=OUT/'local-evidence.json'
    if target.exists():raise FileExistsError(target)
    records={}
    for name in ['symbolic-limit','finite-norm','limit-checks','core-assembly']:
        data=json.loads((OUT/(name+'.json')).read_text())
        assert data['status']=='passed'
        for n,h in data['input_sha256'].items():assert sha(ROOT/n)==h,n
        records[name]=data
    assert records['symbolic-limit']['averaged_coefficients']==['0','0','4*lam**2']
    assert 'gram_over_C_squared' not in records['symbolic-limit']
    finite=records['finite-norm'];limit=records['limit-checks']
    assert finite['general_pair_checks']==1331 and finite['exact_CAR_coordinate_pairs']==31466
    assert len(finite['rows'])==4 and len(limit['records'])==36
    assert sum(r['exact_entries'] for r in limit['finite_checks'])==16
    assert len(records['core-assembly']['gram_B_Linput_O'])==3
    parentpath=ROOT/'results/susskind-memory-central-window/local-evidence.json'
    parent=json.loads(parentpath.read_text());unchanged=[]
    for n,h in parent['hashes'].items():
        if n in ['README.md','PLAN.md','NEXT.md']:continue
        assert sha(ROOT/n)==h,n
        unchanged.append(n)
    log=(OUT/'main.log').read_text(errors='replace')
    assert re.search(r'6\s+pages',log)
    assert not re.search(r'Overfull|undefined|Rerun to get|Fatal error',log)
    pdf=ROOT/'output/pdf/de-sitter-susskind-memory-central-norm.pdf'
    assert sha(pdf)==sha(BUILD/'main.pdf')
    info=subprocess.check_output(['pdfinfo',str(pdf)],text=True,errors='replace')
    assert re.search(r'Pages:\s+6\b',info)
    renders=[BUILD/f'page-{i}.png' for i in range(1,7)]
    assert all(p.exists() for p in renders)
    head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    assert head==parent['current_git_head']
    names=['code/memory_central_norm.py','code/memory_central_norm_asymptotic.py',
           'code/check_memory_central_norm.py','code/check_memory_central_norm_limit.py',
           'code/memory_central_core_assembly.py','code/build_memory_central_norm.ps1',
           'code/record_memory_central_norm.py','manuscript/susskind-memory-central-norm/main.tex',
           'output/pdf/de-sitter-susskind-memory-central-norm.pdf',
           'notes/susskind-memory-central-norm.md','literature/susskind-memory-central-norm-sources.md',
           'journal/2026-11-06-susskind-memory-central-norm.md','README.md','PLAN.md','NEXT.md']
    paths=[ROOT/n for n in names]+sorted(p for p in OUT.rglob('*') if p.is_file())
    result=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),status='Saved locally; no publication',
        scientific_status='General-charge coherent A norm, compact-central uniform input tails, B/L_input/O Gram, strong A-filter replacement and full W2 core expansion. Output-deficit cross, dynamical spectrum extension, charge-tail integration, finite-memory remainder and bulk dictionary remain open.',
        current_git_head=head,git_commit_created=False,git_push_performed=False,archive_created=False,
        parent_manifest_sha256=sha(parentpath),previous_stage_scientific_files_unchanged=unchanged,
        checks=dict(general_rank_pairs=1331,CAR_coordinate_pairs=31466,CAR_systems=4,
                    all_rank_central_cases=36,precisions=[80,120],finite_exact_gram_entries=16),
        environment=dict(python=platform.python_version(),calculations=[d['environment'] for d in records.values()]),
        build=dict(command='powershell -NoProfile -ExecutionPolicy Bypass -File code/build_memory_central_norm.ps1',
                   compiler=log.splitlines()[0],passes_per_build=2,successful_builds=2,final_pages=6,
                   render_command='pdftoppm -r 95 -png output/pdf/de-sitter-susskind-memory-central-norm.pdf tmp/pdfs/susskind-memory-central-norm-build/page',
                   visual_review='All six final page PNGs viewed at 95 dpi. Explicit raw core and symmetric-product definitions added before final double build and full visual review.'),
        failures=[dict(stage='symbolic adaptation',reason='Casimir helper shadowed central charge x; helper renamed.',resolved=True),
                  dict(stage='symbolic polynomial collection',reason='eta=p^2/m is rational in central charge scaling; expanded eta^2 e^2 consistently.',resolved=True),
                  dict(stage='scope cleanup',reason='Removed inherited half-filled four-field Gram from authoritative central result; preserved original source and output.',resolved=True),
                  dict(stage='terminal print',reason='Print referenced removed variable small after scientific checks passed; preserved source/output and reran after print correction.',resolved=True)],
        evidence_scope='Analytic general-charge identities, uniform logarithmic expansion and factorial tails support the operator estimates; finite CAR and high-precision all-rank sums independently verify formulas. Only B/L_input/O Gram is established centrally here, not the S-containing full Gram.',
        process_state='All calculations, compilation and rendering terminal; no active calculation.',
        hashes={p.relative_to(ROOT).as_posix():sha(p) for p in paths},
        render_hashes={p.relative_to(ROOT).as_posix():sha(p) for p in renders},
        build_pass_hashes={f'pass{i}.txt':sha(BUILD/f'pass{i}.txt') for i in [1,2]})
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status='saved locally',pages=6,files=len(paths),
                         previous_scientific_files_unchanged=len(unchanged),git_head=head)))


if __name__=='__main__':main()
