"""Local evidence for mixed memory spectra, with exact inputs and PDF hashes."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,platform,re,subprocess

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-mixed-spectrum'
BUILD=ROOT/'tmp/pdfs/susskind-memory-mixed-spectrum-build'


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    target=OUT/'local-evidence.json'
    if target.exists():raise FileExistsError(target)
    records={}
    for stem in ['cubic-force-cross','charged-lift','matrix-assembly']:
        data=json.loads((OUT/f'{stem}.json').read_text());assert data['status']=='passed'
        for name,digest in data['input_sha256'].items():assert sha(ROOT/name)==digest,name
        records[stem]=data
    cc=records['cubic-force-cross']['certificates']
    assert len(cc)==2 and all(r['templates']==74 and not r['nonzero_templates'] for r in cc)
    assert sum(r['minor_checks'] for r in cc)==5256
    lift=records['charged-lift'];assert lift['covariance_entries']==93076 and lift['direct_identities']==81
    assert sum(r['nonzero'] for d in lift['direct_CAR'] for r in d['rows'])==80
    assert len(records['matrix-assembly']['coefficients'])==9
    parentpath=ROOT/'results/susskind-memory-oriented-kernel/local-evidence.json'
    parent=json.loads(parentpath.read_text());unchanged=[]
    for name,digest in parent['hashes'].items():
        if name in ['README.md','PLAN.md','NEXT.md']:continue
        assert sha(ROOT/name)==digest,name
        unchanged.append(name)
    log=(OUT/'main.log').read_text(errors='replace')
    assert re.search(r'6\s+pages',log)
    assert not re.search(r'Overfull|undefined|Rerun to get|Fatal error',log)
    tex=ROOT/'manuscript/susskind-memory-mixed-spectrum/main.tex'
    assert not any(ord(c)<32 and c not in '\n\r' for c in tex.read_text())
    pdf=ROOT/'output/pdf/de-sitter-susskind-memory-mixed-spectrum.pdf'
    assert sha(pdf)==sha(BUILD/'main.pdf')
    renders=[BUILD/f'page-{i}.png' for i in range(1,7)]
    assert all(p.exists() for p in renders)
    head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    assert head==parent['current_git_head']
    names=['code/memory_mixed_cross.py','code/memory_mixed_lift.py','code/memory_mixed_assembly.py',
           'code/build_memory_mixed_spectrum.ps1','code/record_memory_mixed_spectrum.py',
           'manuscript/susskind-memory-mixed-spectrum/main.tex',
           'output/pdf/de-sitter-susskind-memory-mixed-spectrum.pdf',
           'notes/susskind-memory-mixed-spectrum.md','literature/susskind-memory-mixed-spectrum-sources.md',
           'journal/2026-11-06-susskind-memory-mixed-spectrum.md','README.md','PLAN.md','NEXT.md']
    paths=[ROOT/n for n in names]+sorted(p for p in OUT.rglob('*') if p.is_file())
    data=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),status='Saved locally; no publication',
        scientific_status='Fixed-half-filled R/Y_plus, F/Y_plus and F/Q4 cross measures vanish. The explicitly projected second memory coefficient X is spectrally orthogonal to L and Y; positive three-field matrix and projected polynomial covariance assembled. Raw dynamic descendants, uniform central-charge extension, finite-memory remainder and bulk dictionary remain open.',
        current_git_head=head,git_commit_created=False,git_push_performed=False,archive_created=False,
        parent_manifest_sha256=sha(parentpath),previous_stage_scientific_files_unchanged=unchanged,
        checks=dict(endpoint_templates=148,orbit_words=1752,local_minor_checks=5256,
                    upper_charge_covariance_entries=93076,direct_CAR_identities=81,nonzero_direct_witnesses=80,
                    projected_polynomial_coefficients=9),
        environment=dict(python=platform.python_version(),calculations=[r['environment'] for r in records.values()]),
        build=dict(command='powershell -NoProfile -ExecutionPolicy Bypass -File code/build_memory_mixed_spectrum.ps1',
                   compiler=log.splitlines()[0],passes_per_successful_build=2,successful_builds=2,pages=6,
                   render_command='pdftoppm -r 95 -png output/pdf/de-sitter-susskind-memory-mixed-spectrum.pdf tmp/pdfs/susskind-memory-mixed-spectrum-build/page',
                   visual_review='All six final page PNGs viewed at 95 dpi; equations, references and page boundaries clear. A_minus/A_plus definitions made explicit and lift d distinguished from spectral d_v.'),
        failures=[dict(stage='finite witness coverage',reason='One auxiliary Dbar(gU) is exactly zero for a traceless one-body sample. All identities passed; blanket nonzero assertion was too strong. Zero index recorded; all commutator and quartic witnesses still nonzero.',resolved=True),
                  dict(stage='LaTeX first build',reason='One frac command was serialized as a form-feed. Original source/log saved; corrected.',resolved=True),
                  dict(stage='LaTeX overflow gate',reason='0.57117pt inline determinant overflow; changed to display. Original log saved.',resolved=True)],
        process_state='All calculations, compilation and rendering terminal; no active calculation.',
        evidence_scope='All-order endpoint coverage plus exact lifted covariance and CAR identities; asymptotic cancellation also requires established rank-filter L2 and fixed-word overlap bounds. Three-field result uses an extra projection on the second coefficient and does not establish a raw finite-memory kernel.',
        hashes={p.relative_to(ROOT).as_posix():sha(p) for p in paths},
        render_hashes={p.relative_to(ROOT).as_posix():sha(p) for p in renders},
        build_pass_hashes={f'pass{i}.txt':sha(BUILD/f'pass{i}.txt') for i in [1,2]})
    target.write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(dict(status='saved locally',pages=6,files=len(paths),
                         prior_scientific_files_unchanged=len(unchanged),git_head=head)))


if __name__=='__main__':main()
