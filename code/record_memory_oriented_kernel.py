"""Record the verified local oriented-force and memory spectral stage."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,platform,re,subprocess

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-oriented-kernel'
BUILD=ROOT/'tmp/pdfs/susskind-memory-oriented-kernel-build'


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    target=OUT/'local-evidence.json'
    if target.exists():raise FileExistsError(target)
    records={}
    for stem in ['oriented-marks','force-bridge','reverse-cross','spectral-checks']:
        data=json.loads((OUT/f'{stem}.json').read_text())
        assert data['status']=='passed'
        for name,digest in data['input_sha256'].items():assert sha(ROOT/name)==digest,name
        records[stem]=data
    assert records['oriented-marks']['local_minor_rules']==1688
    assert records['oriented-marks']['finite_comparisons']==2560
    assert records['force-bridge']['finite_identities']==90
    assert records['reverse-cross']['certificate']['templates']==128
    assert records['reverse-cross']['certificate']['words']==6072
    assert not records['reverse-cross']['certificate']['counterexamples']
    assert len(records['spectral-checks']['mixed_moments'])==9
    parentpath=ROOT/'results/susskind-memory-descendants/local-evidence.json'
    parent=json.loads(parentpath.read_text());unchanged=[]
    for name,digest in parent['hashes'].items():
        if name in ['README.md','PLAN.md','NEXT.md']:continue
        assert sha(ROOT/name)==digest,name
        unchanged.append(name)
    log=(OUT/'main.log').read_text(errors='replace')
    assert re.search(r'7\s+pages',log)
    assert not re.search(r'Overfull|undefined|Rerun to get|Fatal error',log)
    tex=ROOT/'manuscript/susskind-memory-oriented-kernel/main.tex'
    assert '\t' not in tex.read_text()
    pdf=ROOT/'output/pdf/de-sitter-susskind-memory-oriented-kernel.pdf'
    assert sha(pdf)==sha(BUILD/'main.pdf')
    renders=[BUILD/f'page-{i}.png' for i in range(1,8)]
    assert all(p.exists() for p in renders)
    head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    assert head==parent['current_git_head']
    names=['code/memory_oriented_marks.py','code/memory_force_orientation_bridge.py',
           'code/memory_oriented_cross.py','code/memory_oriented_spectrum.py',
           'code/build_memory_oriented_kernel.ps1','code/record_memory_oriented_kernel.py',
           'manuscript/susskind-memory-oriented-kernel/main.tex',
           'output/pdf/de-sitter-susskind-memory-oriented-kernel.pdf',
           'notes/susskind-memory-oriented-kernel.md',
           'literature/susskind-memory-oriented-kernel-sources.md',
           'journal/2026-11-06-susskind-memory-oriented-kernel.md','README.md','PLAN.md','NEXT.md']
    paths=[ROOT/n for n in names]+sorted(p for p in OUT.rglob('*') if p.is_file())
    errors=[v['max_absolute_error'] for r in records['spectral-checks']['kernel']['rows']
            for v in r.get('independent_mehler',[])]
    data=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),status='Saved locally; no publication',
        scientific_status='Fixed-half-filled output deficit reduced in L2 to both force orientations. Positive two-energy force matrix and all-order Q4 reverse-force cancellation close the specified eight-direction projection of raw mZ2: V4*mu_q4+V2*mu_q2. Cross with Y_minus and its energy descendants vanishes. Mixed L spectra, central-charge uniformity, finite-memory and bulk dictionary remain open.',
        current_git_head=head,git_commit_created=False,git_push_performed=False,archive_created=False,
        parent_manifest_sha256=sha(parentpath),previous_stage_scientific_files_unchanged=unchanged,
        checks=dict(local_minor_rules=1688,finite_car_comparisons=2560,finite_bridge_identities=90,
                    reverse_endpoint_templates=128,reverse_orbit_words=6072,mixed_low_moments=9,
                    maximum_mehler_jacobi_absolute_difference=max(errors)),
        environment=dict(python=platform.python_version(),calculations=[r['environment'] for r in records.values()]),
        build=dict(command='powershell -NoProfile -ExecutionPolicy Bypass -File code/build_memory_oriented_kernel.ps1',
                   compiler=log.splitlines()[0],passes_per_build=2,successful_builds=2,pages=7,
                   render_command='pdftoppm -r 95 -png output/pdf/de-sitter-susskind-memory-oriented-kernel.pdf tmp/pdfs/susskind-memory-oriented-kernel-build/page',
                   visual_review='All seven pages reviewed individually at 95 dpi. After adding the X/Y_minus spectral cross statement, the final seventh-page layout remains unchanged and final page 5 was reviewed again; all final pages subsequently reviewed before recording.'),
        failures=[],process_state='All calculations, compilation and rendering terminal; no background calculation.',
        evidence_scope='Finite CAR and symbolic certificates verify identities. The L2 asymptotic bridge additionally uses explicit sitewise deficit cancellation and exponential overlap tails. All-order endpoint orbit coverage plus positive matrix measure convergence gives spectra; finite low moments alone do not. Numerical rounding and quadrature are not interval certified.',
        hashes={p.relative_to(ROOT).as_posix():sha(p) for p in paths},
        render_hashes={p.relative_to(ROOT).as_posix():sha(p) for p in renders},
        build_pass_hashes={f'pass{i}.txt':sha(BUILD/f'pass{i}.txt') for i in [1,2]})
    target.write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(dict(status='saved locally',pages=7,files=len(paths),
                         prior_scientific_files_unchanged=len(unchanged),git_head=head)))


if __name__=='__main__':main()
