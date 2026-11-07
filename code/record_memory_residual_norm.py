"""Record fixed-sector norm, remainder proof, and completed PDF evidence."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,platform,re,subprocess

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-residual-norm'
BUILD=ROOT/'tmp/pdfs/susskind-memory-residual-norm-build'


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    target=OUT/'local-evidence.json'
    if target.exists():raise FileExistsError(target)
    records={}
    for stem in ['exchange-casimir','finite-norm','symbolic-limit','limit-checks']:
        data=json.loads((OUT/f'{stem}.json').read_text())
        assert data['status']=='passed',stem
        for name,digest in data['input_sha256'].items():assert sha(ROOT/name)==digest,name
        records[stem]=data
    assert records['exchange-casimir']['ordered_coordinate_pairs']==68544
    assert records['exchange-casimir']['identities']==45
    assert records['finite-norm']['rational_norm_identities']==1200
    assert records['symbolic-limit']['averaged_coefficients']==['0','0','4*lam**2']
    assert records['symbolic-limit']['small_lambda_residual_over_Js_four']=='16'
    assert len(records['limit-checks']['records'])==12
    parentpath=ROOT/'results/susskind-memory-coherent-cross/local-evidence.json'
    parent=json.loads(parentpath.read_text());unchanged=[]
    for name,digest in parent['hashes'].items():
        if name in ['README.md','PLAN.md','NEXT.md']:continue
        assert sha(ROOT/name)==digest,name
        unchanged.append(name)
    log=(OUT/'main.log').read_text(errors='replace')
    assert re.search(r'6\s+pages',log)
    assert not re.search(r'Overfull|undefined|Rerun to get|Fatal error',log)
    pdf=ROOT/'output/pdf/de-sitter-susskind-memory-residual-norm.pdf'
    assert sha(pdf)==sha(BUILD/'main.pdf')
    renders=[BUILD/f'page-{i}.png' for i in range(1,7)]
    assert all(p.exists() for p in renders)
    head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    assert head==parent['current_git_head']
    names=['code/memory_residual_norm.py','code/check_memory_residual_norm.py',
           'code/check_exchange_casimir.py','code/memory_norm_asymptotic.py',
           'code/check_memory_norm_limit.py','code/build_memory_residual_norm.ps1',
           'code/record_memory_residual_norm.py',
           'manuscript/susskind-memory-residual-norm/main.tex',
           'output/pdf/de-sitter-susskind-memory-residual-norm.pdf',
           'notes/susskind-memory-residual-norm.md','literature/susskind-memory-residual-norm-sources.md',
           'journal/2026-11-06-susskind-memory-residual-norm.md','README.md','PLAN.md','NEXT.md']
    paths=[ROOT/n for n in names]+sorted(p for p in OUT.iterdir() if p.is_file())
    data=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),status='Saved locally; no publication',
        scientific_status='Fixed half-filled A norm and the complete four-direction equal-time Gram derived. Output deficit retains a positive residual after projection onto B,A,J_H F. Growing-window proof controls all fixed polynomial input-deficit weights. Output-rank tail, the g_A L2 replacement, central-charge window, raw memory spectrum and bulk dictionary remain open.',
        current_git_head=head,git_commit_created=False,git_push_performed=False,archive_created=False,
        parent_manifest_sha256=sha(parentpath),previous_stage_scientific_files_unchanged=unchanged,
        checks=dict(exchange_casimir_identities=45,ordered_car_coordinate_pairs=68544,
                    lower_charge_contiguous_identities=384,rational_norm_identities=1200,
                    independent_saved_car_gram_entries=27,additional_exact_norms=2,
                    full_rank_high_precision_cases=12),
        environment=dict(python=platform.python_version(),calculations=[r['environment'] for r in records.values()]),
        build=dict(command='powershell -NoProfile -ExecutionPolicy Bypass -File code/build_memory_residual_norm.ps1',
                   compiler=log.splitlines()[0],passes=2,pages=6,
                   render_command='pdftoppm -r 95 -png output/pdf/de-sitter-susskind-memory-residual-norm.pdf tmp/pdfs/susskind-memory-residual-norm-build/page',
                   visual_review='All six final page PNGs viewed individually at 95 dpi; readable equations, references and pagination without clipping or overlap.'),
        failures=[],process_state='All four computations, compilation and rendering terminal; no background calculation.',
        evidence_scope='Finite exact norm and general-charge Casimir lemma, explicit symbolic cancellation, logarithmic-window remainder proof and factorial input tails. Numerical difference/projection use the finite Casimir proxy and retain finite-p deviations. No interval certificate, output-tail theorem or new dynamical particle claim.',
        hashes={p.relative_to(ROOT).as_posix():sha(p) for p in paths},
        render_hashes={p.relative_to(ROOT).as_posix():sha(p) for p in renders},
        build_pass_hashes={f'pass{i}.txt':sha(BUILD/f'pass{i}.txt') for i in [1,2]})
    target.write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(dict(status='saved locally',pages=6,files=len(paths),prior_scientific_files_unchanged=len(unchanged),git_head=head)))


if __name__=='__main__':main()
