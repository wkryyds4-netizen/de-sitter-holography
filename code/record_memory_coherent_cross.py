"""Bind completed coherent cross calculations, proof, build and render evidence."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib, json, platform, re, subprocess

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-coherent-cross'
BUILD=ROOT/'tmp/pdfs/susskind-memory-coherent-cross-build'


def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    target=OUT/'local-evidence.json'
    if target.exists(): raise FileExistsError(target)
    records={}
    for stem in ['transform-checks','charged-map-checks','car-cross','limit-checks']:
        data=json.loads((OUT/f'{stem}.json').read_text())
        assert data['status']=='passed',stem
        for name,digest in data['input_sha256'].items(): assert sha(ROOT/name)==digest,name
        records[stem]=data
    assert records['transform-checks']['embedded_covariance_identities']==1925
    assert records['transform-checks']['signed_unsigned_moments']==556
    assert records['charged-map-checks']['complete_maps']==10
    assert records['charged-map-checks']['complex_entries']==306152
    assert sum(r['coordinate_pairs'] for r in records['car-cross']['rows'])==26227
    assert len(records['limit-checks']['records'])==12
    parentpath=ROOT/'results/susskind-memory-parity-gram/local-evidence.json'
    parent=json.loads(parentpath.read_text());unchanged=[]
    for name,digest in parent['hashes'].items():
        if name in ['README.md','PLAN.md','NEXT.md']:continue
        assert sha(ROOT/name)==digest,name
        unchanged.append(name)
    log=(OUT/'main.log').read_text(errors='replace')
    assert re.search(r'6\s+pages',log)
    assert not re.search(r'Overfull|undefined|Rerun to get|Fatal error',log)
    pdf=ROOT/'output/pdf/de-sitter-susskind-memory-coherent-cross.pdf'
    assert sha(pdf)==sha(BUILD/'main.pdf')
    renders=[BUILD/f'page-{i}.png' for i in range(1,7)]
    assert all(p.exists() for p in renders)
    head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    assert head==parent['current_git_head']
    names=['code/memory_coherent_cross.py','code/check_memory_coherent_transform.py',
           'code/check_memory_coherent_car.py','code/check_memory_coherent_channels.py',
           'code/check_memory_coherent_limit.py','code/build_memory_coherent_cross.ps1',
           'code/record_memory_coherent_cross.py',
           'manuscript/susskind-memory-coherent-cross/main.tex',
           'output/pdf/de-sitter-susskind-memory-coherent-cross.pdf',
           'notes/susskind-memory-coherent-cross.md',
           'literature/susskind-memory-coherent-cross-sources.md',
           'journal/2026-11-06-susskind-memory-coherent-cross.md',
           'README.md','PLAN.md','NEXT.md']
    paths=[ROOT/n for n in names]+sorted(p for p in OUT.iterdir() if p.is_file())
    result=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),
        status='Saved locally; no publication',
        scientific_status='Adjacent-charge coherent transform and exact signed first moments; full same-J half-filled S-A cross derived. Its energy-square-subtracted cross equals the S residual variance. The fixed-sector norm of W2 A minus S, weighted A tails and complete old-span projection remain open.',
        current_git_head=head,git_commit_created=False,git_push_performed=False,archive_created=False,
        parent_manifest_sha256=sha(parentpath),previous_stage_scientific_files_unchanged=unchanged,
        checks=dict(embedded_covariance_identities=1925,signed_unsigned_moments=556,
                    full_charged_maps=10,complex_entries=306152,full_car_coordinate_pairs=26227,
                    contiguous_identities=384,full_rank_high_precision_cases=12),
        environment=dict(python=platform.python_version(),calculations=[r['environment'] for r in records.values()]),
        build=dict(command='powershell -NoProfile -ExecutionPolicy Bypass -File code/build_memory_coherent_cross.ps1',
                   compiler=log.splitlines()[0],passes=2,pages=6,
                   render_command='pdftoppm -r 95 -png output/pdf/de-sitter-susskind-memory-coherent-cross.pdf tmp/pdfs/susskind-memory-coherent-cross-build/page',
                   visual_review='All six final page PNGs viewed individually at 95 dpi. Equations, citations and pagination readable without overlap or clipping.'),
        failures=[dict(record='results/susskind-memory-coherent-cross/verification-issues.txt',
                       scope='Two LaTeX failures: undefined mathscr, then unavailable mathrsfs. Switched to standard mathcal. Scientific calculations all passed first execution.')],
        process_state='All computations, compilation and rendering terminal; CAR session 90090 exited zero.',
        evidence_scope='Analytic finite transform and cross limits with uniform input tails; Casimir proxy converted using a weak positive fixed-sector norm bound. Full-Fock norm is not identified with the exact half-filled norm. Decimal refinement is not an interval certificate; no time-dependent spectrum claimed.',
        reference_hashes={n:sha(ROOT/n) for n in [
            'notes/susskind-singlet-kernel.md','manuscript/susskind-singlet-kernel/main.tex',
            'notes/susskind-memory-parity-gram.md','notes/susskind-memory-two-charge.md']},
        hashes={p.relative_to(ROOT).as_posix():sha(p) for p in paths},
        render_hashes={p.relative_to(ROOT).as_posix():sha(p) for p in renders},
        build_pass_hashes={f'pass{i}.txt':sha(BUILD/f'pass{i}.txt') for i in [1,2]})
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status='saved locally',pages=6,files=len(paths),
                         prior_scientific_files_unchanged=len(unchanged),git_head=head)))


if __name__=='__main__':main()
