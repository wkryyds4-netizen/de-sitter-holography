"""Record completed local descendant and raw-memory equal-time assembly."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,platform,re,subprocess

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-descendants'
BUILD=ROOT/'tmp/pdfs/susskind-memory-descendants-build'


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    target=OUT/'local-evidence.json'
    if target.exists():raise FileExistsError(target)
    records={}
    for stem in ['symbolic-filters','finite-descendants','assembly-certificate']:
        data=json.loads((OUT/f'{stem}.json').read_text())
        assert data['status']=='passed'
        for name,digest in data['input_sha256'].items():assert sha(ROOT/name)==digest,name
        records[stem]=data
    finite=records['finite-descendants']['rows']
    assert sum(r['contraction_identities'] for r in finite)==3600
    assert sum(r['core_raw_identities'] for r in finite)==720
    assert records['assembly-certificate']['fixed_charge_trace_identities']==123
    assert records['assembly-certificate']['rank_mean_identities']==10
    assert len(records['assembly-certificate']['rank_mean_skipped_systems'])==1
    assert len(records['symbolic-filters']['rows'])==5
    parentpath=ROOT/'results/susskind-memory-filter-limit/local-evidence.json'
    parent=json.loads(parentpath.read_text());unchanged=[]
    for name,digest in parent['hashes'].items():
        if name in ['README.md','PLAN.md','NEXT.md']:continue
        assert sha(ROOT/name)==digest,name
        unchanged.append(name)
    log=(OUT/'main.log').read_text(errors='replace')
    assert re.search(r'7\s+pages',log)
    assert not re.search(r'Overfull|undefined|Rerun to get|Fatal error',log)
    tex=ROOT/'manuscript/susskind-memory-descendants/main.tex'
    assert '\t' not in tex.read_text()
    pdf=ROOT/'output/pdf/de-sitter-susskind-memory-descendants.pdf'
    assert sha(pdf)==sha(BUILD/'main.pdf')
    renders=[BUILD/f'page-{i}.png' for i in range(1,8)]
    assert all(p.exists() for p in renders)
    head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    assert head==parent['current_git_head']
    names=['code/memory_descendant_filters.py','code/check_memory_descendants.py',
           'code/check_memory_descendant_assembly.py','code/build_memory_descendants.ps1',
           'code/record_memory_descendants.py','manuscript/susskind-memory-descendants/main.tex',
           'output/pdf/de-sitter-susskind-memory-descendants.pdf',
           'notes/susskind-memory-descendants.md','literature/susskind-memory-descendants-sources.md',
           'journal/2026-11-06-susskind-memory-descendants.md','README.md','PLAN.md','NEXT.md']
    paths=[ROOT/n for n in names]+sorted(p for p in OUT.rglob('*') if p.is_file())
    data=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),status='Saved locally; no publication',
        scientific_status='Raw quartic descendants restored through all ordinary chaos degrees; uniform second-chaos table, marked Wick-block half-filled trace bridge, and the eight-direction equal-time projection of collective raw mZ2 assembled. Residual is Q4 plus the second-chaos output-deficit residual; their energy-resolved cross and complete spectrum remain open.',
        current_git_head=head,git_commit_created=False,git_push_performed=False,archive_created=False,
        parent_manifest_sha256=sha(parentpath),previous_stage_scientific_files_unchanged=unchanged,
        checks=dict(raw_core_identities=720,single_contraction_identities=3600,
                    full_quadratic_tensors=2,larger_half_filled_directions=3,
                    fixed_charge_trace_identities=123,common_rank_scalar_means=10),
        environment=dict(python=platform.python_version(),calculations=[r['environment'] for r in records.values()]),
        build=dict(command='powershell -NoProfile -ExecutionPolicy Bypass -File code/build_memory_descendants.ps1',
                   compiler=log.splitlines()[0],passes_per_build=2,successful_builds=3,pages=7,
                   render_command='pdftoppm -r 95 -png output/pdf/de-sitter-susskind-memory-descendants.pdf tmp/pdfs/susskind-memory-descendants-build/page',
                   visual_review='All seven final page PNGs viewed individually at 95 dpi. A tab-plus-imes typo was detected on an earlier page 5 and corrected to the multiplication symbol before final render review.'),
        failures=[dict(stage='symbolic projection parsing',reason='Saved symbol lambda is a Python keyword; renamed to lam before parsing.',resolved=True),
                  dict(stage='rank scalar mean check',reason='N4p2m3 has a lower-only rank outside the common-rank formula domain; exact CAR checks remain valid. Domain now explicit and skipped case recorded.',resolved=True),
                  dict(stage='visual review',reason='One multiplication command was entered as a tab plus imes; fixed and all final pages reviewed.',resolved=True)],
        process_state='All computations, compilation and rendering terminal; no background calculation.',
        evidence_scope='Finite tensor polarization checks support exact raw contractions; uniform filters and positive tail estimates prove the strong limits. Fixed-charge bridge uses balanced singleton occupations and absolute overlap tails, not numerical extrapolation. Raw mZ2 conclusion is equal-time and for the specified finite span, not a full commutant or spectral theorem.',
        hashes={p.relative_to(ROOT).as_posix():sha(p) for p in paths},
        render_hashes={p.relative_to(ROOT).as_posix():sha(p) for p in renders},
        build_pass_hashes={f'pass{i}.txt':sha(BUILD/f'pass{i}.txt') for i in [1,2]})
    target.write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(dict(status='saved locally',pages=7,files=len(paths),
                         prior_scientific_files_unchanged=len(unchanged),git_head=head)))


if __name__=='__main__':main()
