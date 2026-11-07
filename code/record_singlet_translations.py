"""Save local evidence, preserve earlier science, and verify plotted examples."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,re,subprocess
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-translations'
BUILD=ROOT/'tmp/pdfs/susskind-translations-build'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    target=OUT/'local-evidence.json'
    if target.exists():raise FileExistsError(target)
    previous=[]
    for stage in ['susskind-singlet-moments','susskind-singlet-double-scale','susskind-singlet-kernel','susskind-singlet-multimark','susskind-singlet-otoc','susskind-rindler-packet','susskind-first-order']:
        path=ROOT/f'results/{stage}/local-evidence.json';d=json.loads(path.read_text());kept=[]
        for name,h in d['hashes'].items():
            if name in ['README.md','PLAN.md','NEXT.md']:continue
            assert sha(ROOT/name)==h,name
            kept.append(name)
        previous.append(dict(stage=stage,manifest_sha256=sha(path),unchanged_verified=kept))
    study=json.loads((OUT/'study.json').read_text());assert study['status']=='passed'
    plot=json.loads((OUT/'plot-data.json').read_text())
    for d in [study,plot]:
        for p,h in d['input_sha256'].items():assert sha(ROOT/p)==h,p
    idx=plot['a'].index(.5);plot_checks=[]
    for ref,curve in zip(study['mass_comparison'],plot['curves']):
        assert float(ref['m'])==curve['m']
        z=complex(ref['spacetime_value'].replace(' ',''));v=complex(curve['real'][idx],curve['imag'][idx])
        err=abs(v-z);assert err<2e-14
        plot_checks.append(dict(m=curve['m'],absolute_difference=err))
    assert len(study['symbolic']['operator_checks'])==9 and len(study['kernels'])==8
    log=(OUT/'main.log').read_text(errors='replace')
    assert not re.search(r'Overfull|undefined|Rerun to get|Fatal error',log)
    assert '(6 pages,' in log
    assert sha(ROOT/'output/pdf/de-sitter-susskind-translations.pdf')==sha(BUILD/'main.pdf')
    renders=[BUILD/f'final-page-{i}.png' for i in range(1,7)]
    assert all(p.exists() for p in renders)
    for i in [1,2,3,4,6]:assert sha(BUILD/f'page-{i}.png')==sha(renders[i-1])
    names=['code/singlet_translations.py','code/check_singlet_translations.py','code/plot_singlet_translations.py',
        'code/build_singlet_translations.ps1','code/record_singlet_translations.py',
        'manuscript/susskind-translations/main.tex','output/pdf/de-sitter-susskind-translations.pdf',
        'notes/susskind-translations.md','literature/susskind-translations-sources.md','journal/2026-11-06-susskind-translations.md',
        'README.md','PLAN.md','NEXT.md']
    files=[ROOT/p for p in names]+sorted(p for p in OUT.iterdir() if p.is_file())
    data=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),status='Local stage saved; no publication or global freeze',
        scientific_status='Explicit candidate positive-energy extensions, spacetime kernels and finite-q seed mismatch; microscopic translations and full de Sitter duality remain open.',
        previous_stages=previous,current_git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        git_commit_created=False,git_push_performed=False,archive_created=False,global_project_audit_performed=False,
        plot_checks=plot_checks,checks=dict(exact_differential_identities=9,quadrature_cases=8,boost_cases=4),
        build=dict(command='powershell -NoProfile -ExecutionPolicy Bypass -File code/build_singlet_translations.ps1',compiler=log.splitlines()[0],passes=2,pages=6,
            visual_review='All six 1200px pages inspected; final pages 1/2/3/4/6 match inspected images, final page 5 re-inspected after trace argument addition.'),
        sources='Primary papers read online; precise scope and URLs in literature note; no new source-PDF cache.',
        process_state='All numerical/build/render sessions terminal',
        hashes={p.relative_to(ROOT).as_posix():sha(p) for p in files},local_render_hashes={p.relative_to(ROOT).as_posix():sha(p) for p in renders},
        build_pass_hashes={f'pass{i}.txt':sha(BUILD/f'pass{i}.txt') for i in [1,2]})
    target.write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(dict(status='saved locally',files=len(files),previous_scientific_files=[len(d['unchanged_verified']) for d in previous],pages=6,git_head=data['current_git_head'],plot_checks=plot_checks)))
if __name__=='__main__':main()
