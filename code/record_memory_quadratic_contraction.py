"""Record the local contraction stage; preserve prior scientific evidence."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,re,subprocess,platform

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-quadratic-contraction'
BUILD=ROOT/'tmp/pdfs/susskind-memory-quadratic-contraction-build'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    target=OUT/'local-evidence.json'
    if target.exists():raise FileExistsError(target)
    parent=json.loads((ROOT/'results/susskind-memory-quartic/local-evidence.json').read_text())
    previous=[]
    for stage in [r['stage'] for r in parent['previous_stages']]+['susskind-memory-quartic']:
        path=ROOT/f'results/{stage}/local-evidence.json';old=json.loads(path.read_text());kept=[]
        for name,digest in old['hashes'].items():
            if name in ['README.md','PLAN.md','NEXT.md']:continue
            assert sha(ROOT/name)==digest,name
            kept.append(name)
        previous.append(dict(stage=stage,manifest_sha256=sha(path),unchanged_verified=kept))
    records={}
    for name in ['finite-contractions','scalar-removed-checks','onebody-controls','linear-input-audit']:
        r=json.loads((OUT/f'{name}.json').read_text());assert r['status']=='passed'
        for path,digest in r['input_sha256'].items():assert sha(ROOT/path)==digest,path
        records[name]=r
    a,b,c,d=[records[n] for n in ['finite-contractions','scalar-removed-checks','onebody-controls','linear-input-audit']]
    assert sum(r['quadratic_coefficients'] for r in a['rows'])==5761
    assert sum(r['independent_polarization_directions'] for r in b['rows'])==807
    assert [r['residual_norm_squared'] for r in b['rows']]==['0','0','248945355/14903']
    assert b['rows'][2]['residual_commutator_norm_squared']=='12991771132439678939/2220994090'
    assert c['identities']==56 and sum(r['directions'] for r in d['rows'])==148
    log=(OUT/'main.log').read_text(errors='replace')
    assert re.search(r'5\s+pages',log)
    assert not re.search(r'Overfull|undefined|Rerun to get|Fatal error',log)
    pdf=ROOT/'output/pdf/de-sitter-susskind-memory-quadratic-contraction.pdf'
    assert sha(pdf)==sha(BUILD/'main.pdf')
    renders=[BUILD/f'page-{i}.png' for i in range(1,6)]
    assert all(p.exists() for p in renders)
    head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    assert head==parent['current_git_head']
    names=['code/memory_quadratic_contraction.py','code/check_quadratic_contraction.py',
           'code/check_quadratic_onebody.py','code/audit_quadratic_inputs.py',
           'code/build_memory_quadratic_contraction.ps1','code/record_memory_quadratic_contraction.py',
           'manuscript/susskind-memory-quadratic-contraction/main.tex',
           'output/pdf/de-sitter-susskind-memory-quadratic-contraction.pdf',
           'notes/susskind-memory-quadratic-contraction.md',
           'literature/susskind-memory-quadratic-contraction-sources.md',
           'journal/2026-11-06-susskind-memory-quadratic-contraction.md','README.md','PLAN.md','NEXT.md']
    files=[ROOT/n for n in names]+sorted(p for p in OUT.iterdir() if p.is_file())
    result=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),status='Saved locally; no publication or global freeze',
        scientific_status='Exact finite-N W2 memory-core contraction map, energy-quadratic projection and sample-scalar removal; explicit finite dynamic residual beyond old directions. DS residual strength, raw Z2 spectrum and finite-memory/bulk dictionary remain open.',
        previous_stages=previous,current_git_head=head,git_commit_created=False,git_push_performed=False,
        archive_created=False,global_project_audit_performed=False,
        checks=dict(quadratic_coefficient_blocks=5761,independent_quartic_polarization=807,
                    onebody_identities=56,independent_cubic_input_directions=148),
        environment=dict(python=platform.python_version(),calculation_environments=[r['environment'] for r in records.values()]),
        build=dict(command='powershell -NoProfile -ExecutionPolicy Bypass -File code/build_memory_quadratic_contraction.ps1',
                   compiler=log.splitlines()[0],passes=2,pages=5,
                   visual_review='All five final pages viewed at 95 dpi; no overflow, clipping, overlap, or missing citations.',
                   render_command='pdftoppm -r 95 -png output/pdf/de-sitter-susskind-memory-quadratic-contraction.pdf tmp/pdfs/susskind-memory-quadratic-contraction-build/page'),
        retained_issues=['Initial two-pass PDF had a 2.48915pt paragraph overflow. Source/log retained; sentence split and final two-pass build clean.',
                         'Positive finite Gram residuals can be random scalar operators: sample-wise trace removal closes the two smaller examples. Only the N5p2m3 residual here has a positive explicit commutator witness.'],
        evidence_scope='Finite-sector exact coordinate contractions and rational Gram, checked by independent CAR polarization. No Monte Carlo, no DS extrapolation, no claim of complete commutant subtraction.',
        process_state='All calculations, compilation and rendering terminal; scalar checker session46093 exited zero.',
        hashes={p.relative_to(ROOT).as_posix():sha(p) for p in files},
        local_render_hashes={p.relative_to(ROOT).as_posix():sha(p) for p in renders},
        build_pass_hashes={f'pass{i}.txt':sha(BUILD/f'pass{i}.txt') for i in [1,2]})
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status='saved locally',files=len(files),pages=5,
                         prior_scientific_files=sum(len(r['unchanged_verified']) for r in previous),git_head=head)))


if __name__=='__main__':main()
