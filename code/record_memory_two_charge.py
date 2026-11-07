"""Bind exact filters, asymptotics, proof note, and rendering to a local manifest."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,re,subprocess,platform
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-two-charge'
BUILD=ROOT/'tmp/pdfs/susskind-memory-two-charge-build'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    target=OUT/'local-evidence.json'
    if target.exists():raise FileExistsError(target)
    records={}
    for stem in ['filter-checks','transfer-checks','symbolic-filter','limit-checks','uniform-prefactors']:
        data=json.loads((OUT/f'{stem}.json').read_text());assert data['status']=='passed'
        for name,digest in data['input_sha256'].items():assert sha(ROOT/name)==digest,name
        records[stem]=data
    assert records['filter-checks']['identities']==522
    tr=records['transfer-checks']['rows']
    assert sum(r['full_rectangular_map_identities'] for r in tr)==87
    assert sum(r['complex_entry_equalities'] for r in tr)==5470368
    assert len(records['limit-checks']['records'])==24
    parentpath=ROOT/'results/susskind-memory-charge-lift/local-evidence.json'
    parent=json.loads(parentpath.read_text());unchanged=[]
    for name,digest in parent['hashes'].items():
        if name in ['README.md','PLAN.md','NEXT.md']:continue
        assert sha(ROOT/name)==digest,name
        unchanged.append(name)
    log=(OUT/'main.log').read_text(errors='replace')
    assert re.search(r'6\s+pages',log)
    assert not re.search(r'Overfull|undefined|Rerun to get|Fatal error',log)
    pdf=ROOT/'output/pdf/de-sitter-susskind-memory-two-charge.pdf'
    assert sha(pdf)==sha(BUILD/'main.pdf')
    renders=[BUILD/f'page-{i}.png' for i in range(1,7)]
    assert all(p.exists() for p in renders)
    head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    assert head==parent['current_git_head']
    names=['code/memory_two_charge_filter.py','code/check_memory_two_charge_filter.py',
        'code/check_charged_covariance_transfer.py','code/memory_two_charge_asymptotic.py',
        'code/check_memory_filter_limit.py','code/check_memory_filter_uniform.py',
        'code/build_memory_two_charge.ps1','code/record_memory_two_charge.py',
        'manuscript/susskind-memory-two-charge/main.tex',
        'output/pdf/de-sitter-susskind-memory-two-charge.pdf',
        'notes/susskind-memory-two-charge.md','literature/susskind-memory-two-charge-sources.md',
        'journal/2026-11-06-susskind-memory-two-charge.md','README.md','PLAN.md','NEXT.md']
    paths=[ROOT/n for n in names]+sorted(p for p in OUT.iterdir() if p.is_file())
    result=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),status='Saved locally; no publication',
        scientific_status='General finite-N two-charge core filters proved; half-filled H-squared part admits an L2 linear-deficit approximation. Weighted A tail, parity-weighted output-deficit cross spectrum, full charge window, raw Z2 and bulk dictionary remain open.',
        current_git_head=head,git_commit_created=False,git_push_performed=False,archive_created=False,
        parent_manifest_sha256=sha(parentpath),previous_stage_scientific_files_unchanged=unchanged,
        checks=dict(rank_contraction_identities=522,full_linear_maps=87,complex_entry_equalities=5470368,
            high_precision_cases=24,symbolic_projection_identities=3),
        environment=dict(python=platform.python_version(),calculations=[r['environment'] for r in records.values()]),
        build=dict(command='powershell -NoProfile -ExecutionPolicy Bypass -File code/build_memory_two_charge.ps1',
            compiler=log.splitlines()[0],passes=2,pages=6,
            render_command='pdftoppm -r 95 -png output/pdf/de-sitter-susskind-memory-two-charge.pdf tmp/pdfs/susskind-memory-two-charge-build/page',
            visual_review='All six final page PNGs viewed after the final rebuild at 95 dpi. Equations, references, and pagination readable without clipping or overlapping.'),
        failures=[],process_state='All computations, compilation, and rendering terminal.',
        evidence_scope='Analytic transfer proof plus exhaustive finite linear-map checks; rank-product checks are selected pure-rank realizations. Fixed-deficit asymptotics and a uniform-tail proof concern the H-squared contribution only. Decimal refinements are not interval certificates.',
        hashes={p.relative_to(ROOT).as_posix():sha(p) for p in paths},
        render_hashes={p.relative_to(ROOT).as_posix():sha(p) for p in renders},
        build_pass_hashes={f'pass{i}.txt':sha(BUILD/f'pass{i}.txt') for i in [1,2]})
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status='saved locally',files=len(paths),pages=6,prior_scientific_files_unchanged=len(unchanged),git_head=head)))

if __name__=='__main__':main()
