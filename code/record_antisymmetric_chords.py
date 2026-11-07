"""Record this local stage after all six rendered pages have been inspected."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,re,subprocess
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-antisymmetric'
BUILD=ROOT/'tmp/pdfs/susskind-antisymmetric-build'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    target=OUT/'local-evidence.json'
    if target.exists():raise FileExistsError(target)
    previous=[]
    for stage in ['susskind-singlet-moments','susskind-singlet-double-scale','susskind-singlet-kernel','susskind-singlet-multimark','susskind-singlet-otoc','susskind-rindler-packet','susskind-first-order','susskind-translations']:
        path=ROOT/f'results/{stage}/local-evidence.json';d=json.loads(path.read_text());kept=[]
        for name,h in d['hashes'].items():
            if name in ['README.md','PLAN.md','NEXT.md']:continue
            assert sha(ROOT/name)==h,name
            kept.append(name)
        previous.append(dict(stage=stage,manifest_sha256=sha(path),unchanged_verified=kept))
    study=json.loads((OUT/'validation.json').read_text());assert study['status']=='passed'
    for p,h in study['input_sha256'].items():assert sha(ROOT/p)==h,p
    assert len(study['cases'])==6 and len(study['symbolic_laguerre'])==9
    log=(OUT/'main.log').read_text(errors='replace')
    assert not re.search(r'Overfull|undefined|Rerun to get|Fatal error',log)
    assert '(6 pages,' in log
    assert sha(ROOT/'output/pdf/de-sitter-susskind-antisymmetric.pdf')==sha(BUILD/'main.pdf')
    renders=[BUILD/f'page-{i}.png' for i in range(1,7)]
    assert all(p.exists() for p in renders)
    names=['code/antisymmetric_chord_probe.py','code/check_antisymmetric_chords.py','code/build_antisymmetric_chords.ps1','code/record_antisymmetric_chords.py',
        'manuscript/susskind-antisymmetric/main.tex','output/pdf/de-sitter-susskind-antisymmetric.pdf',
        'notes/susskind-antisymmetric.md','literature/susskind-antisymmetric-sources.md','journal/2026-11-06-susskind-antisymmetric.md',
        'README.md','PLAN.md','NEXT.md']
    files=[ROOT/p for p in names]+sorted(p for p in OUT.iterdir() if p.is_file())
    head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    assert head=='714d038bb605cfcffdc351e4dcc62aa95237d25f'
    data=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),status='Local stage saved; no publication or global freeze',
        scientific_status='Fixed-degree chord contraction and finite-h positive Q_h, with seed spectral convergence. Multiparticle dynamics and de Sitter duality remain open.',
        previous_stages=previous,current_git_head=head,git_commit_created=False,git_push_performed=False,archive_created=False,global_project_audit_performed=False,
        checks=dict(symbolic_laguerre_identities=18,parameter_cases=6,decimal_precisions=[80,110],exact_probe_degrees=list(range(1,7))),
        build=dict(command='powershell -NoProfile -ExecutionPolicy Bypass -File code/build_antisymmetric_chords.ps1',compiler=log.splitlines()[0],passes=2,pages=6,visual_review='All six final 1200px page images inspected; no clipped formulas or references.'),
        process_state='All calculation/build/render sessions terminal',sources='Primary papers read online; scope in literature note.',
        hashes={p.relative_to(ROOT).as_posix():sha(p) for p in files},
        local_render_hashes={p.relative_to(ROOT).as_posix():sha(p) for p in renders},
        build_pass_hashes={f'pass{i}.txt':sha(BUILD/f'pass{i}.txt') for i in [1,2]})
    target.write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(dict(status='saved locally',files=len(files),prior_stage_counts=[len(d['unchanged_verified']) for d in previous],pages=6,git_head=head)))
if __name__=='__main__':main()
