"""Save the two-matter local stage and verify preceding scientific hashes."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,re,subprocess
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-two-matter'
BUILD=ROOT/'tmp/pdfs/susskind-two-matter-build'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    target=OUT/'local-evidence.json'
    if target.exists():raise FileExistsError(target)
    previous=[]
    for stage in ['susskind-singlet-moments','susskind-singlet-double-scale','susskind-singlet-kernel','susskind-singlet-multimark','susskind-singlet-otoc','susskind-rindler-packet','susskind-first-order','susskind-translations','susskind-antisymmetric']:
        path=ROOT/f'results/{stage}/local-evidence.json';d=json.loads(path.read_text());kept=[]
        for name,h in d['hashes'].items():
            if name in ['README.md','PLAN.md','NEXT.md']:continue
            assert sha(ROOT/name)==h,name
            kept.append(name)
        previous.append(dict(stage=stage,manifest_sha256=sha(path),unchanged_verified=kept))
    study=json.loads((OUT/'validation.json').read_text());assert study['status']=='passed'
    for p,h in study['input_sha256'].items():assert sha(ROOT/p)==h,p
    tensor=json.loads((OUT/'tensor-check.json').read_text());assert tensor['status']=='passed'
    assert sha(ROOT/'code/check_two_matter_tensor.py')==tensor['input_sha256']
    assert len(study['symbolic']['checks'])==38 and len(tensor['checks'])==3
    log=(OUT/'main.log').read_text(errors='replace')
    assert not re.search(r'Overfull|undefined|Rerun to get|Fatal error',log)
    assert '(6 pages,' in log
    assert sha(ROOT/'output/pdf/de-sitter-susskind-two-matter.pdf')==sha(BUILD/'main.pdf')
    renders=[BUILD/f'final-page-{i}.png' for i in range(1,7)]
    assert all(p.exists() for p in renders)
    for n in [1,2,5,6]:assert sha(BUILD/f'page-{n}.png')==sha(BUILD/f'final-page-{n}.png')
    names=['code/two_matter_chords.py','code/check_two_matter_chords.py','code/check_two_matter_tensor.py','code/build_two_matter_chords.ps1','code/record_two_matter_chords.py',
        'manuscript/susskind-two-matter/main.tex','output/pdf/de-sitter-susskind-two-matter.pdf',
        'notes/susskind-two-matter.md','literature/susskind-two-matter-sources.md','journal/2026-11-06-susskind-two-matter.md',
        'README.md','PLAN.md','NEXT.md']
    files=[ROOT/p for p in names]+sorted(p for p in OUT.iterdir() if p.is_file())
    head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    assert head=='714d038bb605cfcffdc351e4dcc62aa95237d25f'
    data=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),status='Local stage saved; no publication or global freeze',
        scientific_status='Exact first relative-channel commutator decomposition and ordered two-matter positive generator; independent single-meson spectrum and full de Sitter duality remain open.',
        previous_stages=previous,current_git_head=head,git_commit_created=False,git_push_performed=False,archive_created=False,global_project_audit_performed=False,
        checks=dict(exact_identities=41,parameter_cases=3,decimal_precisions=[85,115],primary_degrees=list(range(5))),
        build=dict(command='powershell -NoProfile -ExecutionPolicy Bypass -File code/build_two_matter_chords.ps1',compiler=log.splitlines()[0],passes=2,pages=6,visual_review='All six pages inspected; final pages 3/4 re-inspected after notation corrections and remaining pages pixel-identical.'),
        process_state='All numerical/build/render sessions terminal',sources='Primary sources read online; exact scope in literature note.',
        hashes={p.relative_to(ROOT).as_posix():sha(p) for p in files},local_render_hashes={p.relative_to(ROOT).as_posix():sha(p) for p in renders},
        build_pass_hashes={f'pass{i}.txt':sha(BUILD/f'pass{i}.txt') for i in [1,2]})
    target.write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(dict(status='saved locally',files=len(files),prior_stage_counts=[len(d['unchanged_verified']) for d in previous],pages=6,git_head=head)))
if __name__=='__main__':main()
