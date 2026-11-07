"""Save the charge-lift result without committing, pushing, or rewriting old evidence."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,re,subprocess,platform

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-charge-lift'
BUILD=ROOT/'tmp/pdfs/susskind-memory-charge-lift-build'
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    target=OUT/'local-evidence.json'
    if target.exists():raise FileExistsError(target)
    data=json.loads((OUT/'final-checks.json').read_text());assert data['status']=='passed'
    for name,digest in data['input_sha256'].items():assert sha(ROOT/name)==digest,name
    preliminary=json.loads((OUT/'checks.json').read_text())
    for name,digest in preliminary['input_sha256'].items():
        path=OUT/'preliminary-checker.py' if name=='code/check_memory_charge_lift.py' else ROOT/name
        assert sha(path)==digest,name
    parentpath=ROOT/'results/susskind-memory-rank-blocks/local-evidence.json'
    parent=json.loads(parentpath.read_text());unchanged=[]
    for name,digest in parent['hashes'].items():
        if name in ['README.md','PLAN.md','NEXT.md']:continue
        assert sha(ROOT/name)==digest,name
        unchanged.append(name)
    assert len(data['family_cases'])==15 and len(data['charge_lift_checks'])==8
    assert data['witness']['P1_Q_core']=='99D'
    assert data['witness']['independent_full_Gaussian_polarization_coordinates']==225
    log=(OUT/'main.log').read_text(errors='replace')
    assert re.search(r'4\s+pages',log)
    assert not re.search(r'Overfull|undefined|Rerun to get|Fatal error',log)
    pdf=ROOT/'output/pdf/de-sitter-susskind-memory-charge-lift.pdf'
    assert sha(pdf)==sha(BUILD/'main.pdf')
    renders=[BUILD/f'page-{i}.png' for i in range(1,5)]
    assert all(p.exists() for p in renders)
    head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    assert head==parent['current_git_head']
    names=['code/check_memory_charge_lift.py','code/probe_memory_half_filling.py',
        'code/build_memory_charge_lift.ps1','code/record_memory_charge_lift.py',
        'manuscript/susskind-memory-charge-lift/main.tex',
        'output/pdf/de-sitter-susskind-memory-charge-lift.pdf',
        'notes/susskind-memory-charge-lift.md','literature/susskind-memory-charge-lift-sources.md',
        'journal/2026-11-06-susskind-memory-charge-lift.md','README.md','PROJECT.md','PLAN.md','NEXT.md']
    paths=[ROOT/n for n in names]+sorted(p for p in OUT.iterdir() if p.is_file())
    result=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),
        status='Saved locally; no publication',
        scientific_status='General scalar multiplication filter disproved; coherent adjacent-charge reduced amplitudes required. Full Q amplitudes, annealed DS residual, raw Z2 and bulk dictionary remain open.',
        current_git_head=head,git_commit_created=False,git_push_performed=False,archive_created=False,
        parent_manifest_sha256=sha(parentpath),previous_stage_scientific_files_unchanged=unchanged,
        checks=dict(analytic_family_systems=15,charge_lift_cases=8,full_rank_parity_projectors=4,
            independent_core_polarization_coordinates=225),
        environment=dict(python=platform.python_version(),calculation=data['environment']),
        build=dict(command='powershell -NoProfile -ExecutionPolicy Bypass -File code/build_memory_charge_lift.ps1',
            compiler=log.splitlines()[0],passes=2,pages=4,
            render_command='pdftoppm -r 95 -png output/pdf/de-sitter-susskind-memory-charge-lift.pdf tmp/pdfs/susskind-memory-charge-lift-build/page',
            visual_review='All four final page PNGs viewed at 95 dpi; equations, references, and pagination readable without clipping or overlaps.'),
        retained_issues=['Ad hoc python -c import failed from root; rerun from code directory succeeded.',
            'Initial passing checks retained with exact source snapshot; final checker adds analytic norm-count validation. No scientific tolerance changes.'],
        process_state='All calculations, compilation, and rendering terminal; no active process.',
        evidence_scope='Finite exact identities and an analytic counterexample family. Does not establish general block multiplicity or annealed double-scaled memory norm.',
        hashes={p.relative_to(ROOT).as_posix():sha(p) for p in paths},
        render_hashes={p.relative_to(ROOT).as_posix():sha(p) for p in renders},
        build_pass_hashes={f'pass{i}.txt':sha(BUILD/f'pass{i}.txt') for i in [1,2]})
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status='saved locally',files=len(paths),pages=4,
        prior_scientific_files_unchanged=len(unchanged),git_head=head)))

if __name__=='__main__':main()
