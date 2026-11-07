"""Record the charged bridge and completed subtracted-cubic spectrum locally."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,re,subprocess,platform
import numpy,mpmath,sympy,matplotlib
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-filtered-spectrum'
BUILD=ROOT/'tmp/pdfs/susskind-filtered-spectrum-build'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    target=OUT/'local-evidence.json'
    if target.exists():raise FileExistsError(target)
    parent=json.loads((ROOT/'results/susskind-cubic-cross/local-evidence.json').read_text())
    stages=[r['stage'] for r in parent['previous_stages']]+['susskind-cubic-cross'];previous=[]
    for stage in stages:
        path=ROOT/f'results/{stage}/local-evidence.json';record=json.loads(path.read_text());kept=[]
        for name,h in record['hashes'].items():
            if name in ['README.md','PLAN.md','NEXT.md']:continue
            assert sha(ROOT/name)==h,name;kept.append(name)
        previous.append(dict(stage=stage,manifest_sha256=sha(path),unchanged_verified=kept))
    data={}
    for name in ['charged-bridge','spectrum','spectral-checks']:
        record=json.loads((OUT/f'{name}.json').read_text());data[name]=record;assert record['status']=='passed'
        for path,h in record['input_sha256'].items():assert sha(ROOT/path)==h,path
    bridge=data['charged-bridge'];spectrum=data['spectrum'];spectral=data['spectral-checks']
    assert len(bridge['finite'])==3 and len(bridge['rank_records'])==12
    assert sum(len(r['covariance_blocks']) for r in bridge['finite'])==8
    assert sum(c['entries'] for r in bridge['finite'] for c in r['covariance_blocks'])==288772
    assert sum(r['marked_product_checks'] for r in bridge['finite'])==15
    assert len(spectrum['records'])==5 and len(spectrum['independent_mehler'])==2
    assert len(spectral['chord_moments'])==9 and len(spectral['fourier'])==4
    assert spectral['exact_leading_norm']=='176' and spectral['exact_leading_second_moment']=='8832'
    failure=json.loads((OUT/'initial-check-failure.json').read_text())
    assert failure['status']=='failed' and (ROOT/failure['original_code']).exists()
    assert 'Overfull' in (OUT/'build-initial.log').read_text(errors='replace')
    log=(OUT/'main.log').read_text(errors='replace')
    assert not re.search(r'Overfull|undefined|Rerun to get|Fatal error',log) and '(7 pages,' in log
    pdf=ROOT/'output/pdf/de-sitter-susskind-filtered-spectrum.pdf'
    assert sha(pdf)==sha(BUILD/'main.pdf')
    renders=[BUILD/f'final-page-{i}.png' for i in range(1,8)];assert all(p.exists() for p in renders)
    for i in range(1,7):assert sha(BUILD/f'page-{i}.png')==sha(BUILD/f'final-page-{i}.png')
    names=['code/singlet_filtered_bridge.py','code/run_filtered_spectrum.py','code/check_filtered_spectrum.py',
        'code/build_filtered_spectrum.ps1','code/record_filtered_spectrum.py',
        'manuscript/susskind-filtered-spectrum/main.tex','output/pdf/de-sitter-susskind-filtered-spectrum.pdf',
        'notes/susskind-filtered-spectrum.md','literature/susskind-filtered-spectrum-sources.md',
        'journal/2026-11-06-susskind-filtered-spectrum.md','README.md','PLAN.md','NEXT.md']
    files=[ROOT/n for n in names]+sorted(p for p in OUT.iterdir() if p.is_file())
    head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    assert head=='714d038bb605cfcffdc351e4dcc62aa95237d25f'
    record=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),status='Local stage saved; no publication or global freeze',
        scientific_status='Charged covariance lift and marked-product identity prove the remaining F-Y/F-R two-energy cross measures vanish in DS. Complete positive measure for U=M3-PN and its fixed-time small-lambda limit. Mixed higher-point interactions and bulk dictionary remain open.',
        previous_stages=previous,current_git_head=head,git_commit_created=False,git_push_performed=False,
        archive_created=False,global_project_audit_performed=False,
        checks=dict(charged_blocks=8,charged_entries=288772,marked_products=15,full_rank_cases=12,
            spectrum_couplings=5,independent_mehler_couplings=2,exact_chord_moments=9,fourier_points=4),
        environment=dict(python=platform.python_version(),numpy=numpy.__version__,mpmath=mpmath.__version__,
            sympy=sympy.__version__,matplotlib=matplotlib.__version__),
        build=dict(command='powershell -NoProfile -ExecutionPolicy Bypass -File code/build_filtered_spectrum.ps1',
            compiler=log.splitlines()[0],passes=2,pages=7,
            visual_review='All seven pages inspected; final page 7 reinspected after additional validation summary, pages 1-6 pixel-identical.'),
        error_scope='Analytic weighted-Jacobi Taylor truncation; full-rank, Fourier and Mehler precision empirical, no interval rounding certification',
        retained_failures=['Initial 1.11525pt typesetting overflow, fixed by shortening a sentence',
            'Structural equality in new spectral checker, fixed by exact expansion of residual; original code retained'],
        process_state='All computation, compilation and rendering processes terminal',
        hashes={p.relative_to(ROOT).as_posix():sha(p) for p in files},
        local_render_hashes={p.relative_to(ROOT).as_posix():sha(p) for p in renders},
        build_pass_hashes={f'pass{i}.txt':sha(BUILD/f'pass{i}.txt') for i in [1,2]})
    target.write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(dict(status='saved locally',files=len(files),
        prior_scientific_files=sum(len(r['unchanged_verified']) for r in previous),pages=7,git_head=head)))
if __name__=='__main__':main()
