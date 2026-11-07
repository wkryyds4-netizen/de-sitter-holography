"""Record this local research stage and verify preceding scientific files."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,re,subprocess

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-singlet-kernel'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    target=OUT/'local-evidence.json'
    if target.exists():raise FileExistsError(target)
    prior=[]
    for stage in ['susskind-singlet-moments','susskind-singlet-double-scale']:
        record=ROOT/f'results/{stage}/local-evidence.json'
        data=json.loads(record.read_text());unchanged=[];updated=[]
        for p,digest in data['hashes'].items():
            if p in ['README.md','PLAN.md','NEXT.md']:
                updated.append(p)
            else:
                assert sha(ROOT/p)==digest,p
                unchanged.append(p)
        prior.append(dict(stage=stage,record_sha256=sha(record),unchanged_verified=unchanged,
                          intentionally_updated_navigation=updated))
    for name in ['exact-checks.json','nonzero-words.json','kernel-numerics.json']:
        d=json.loads((OUT/name).read_text())
        assert d['status']=='passed'
        for p,h in d['input_sha256'].items():assert sha(ROOT/p)==h,p
    exact=json.loads((OUT/'exact-checks.json').read_text())
    extra=json.loads((OUT/'nonzero-words.json').read_text())
    assert exact['car_checks']['scalar_comparisons']+extra['scalar_comparisons']==4320
    assert extra['nonzero']==262 and extra['negative']==23
    build=ROOT/'tmp/pdfs/susskind-singlet-kernel-build'
    log=(OUT/'main.log').read_text(errors='replace')
    assert not re.search(r'Overfull|undefined|Fatal error|Rerun to get',log)
    assert '(8 pages,' in log
    for i in [1,2]:
        raw=(build/f'pass{i}.txt').read_bytes()
        enc='utf-16' if raw.startswith((b'\xff\xfe',b'\xfe\xff')) else 'utf-8'
        assert '(8 pages,' in raw.decode(enc,errors='replace')
    assert sha(build/'main.pdf')==sha(ROOT/'output/pdf/de-sitter-susskind-singlet-kernel.pdf')
    renders=[build/f'page-{i}.png' for i in range(1,9)]
    assert all(p.exists() for p in renders)
    files=[ROOT/p for p in [
        'code/singlet_spectral_graphs.py','code/check_singlet_spectral_kernel.py',
        'code/check_singlet_kernel_nonzero.py','code/run_singlet_kernel.py',
        'code/build_singlet_kernel.ps1','code/record_singlet_kernel_local.py',
        'manuscript/susskind-singlet-kernel/main.tex','output/pdf/de-sitter-susskind-singlet-kernel.pdf',
        'notes/susskind-singlet-kernel.md','literature/susskind-singlet-kernel-sources.md',
        'journal/2026-11-06-susskind-singlet-kernel.md','README.md','PLAN.md','NEXT.md']]
    files+=sorted(p for p in OUT.iterdir() if p.is_file())
    payload=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),
        status='local stage recorded; not frozen or published',
        scientific_status='All fixed mixed-moment endpoint exchange proof gives the specified singlet two-point spectral measure at fixed-positive-lambda double scaling. Bulk mass dictionary and higher correlators remain open.',
        current_git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        git_commit_created=False,git_push_performed=False,archive_created=False,
        global_project_audit_performed=False,previous_stages=prior,
        build=dict(command='powershell -NoProfile -ExecutionPolicy Bypass -File code/build_singlet_kernel.ps1',
                   compiler=log.splitlines()[0],passes=2,pages=8,
                   visual_review='All eight final 1200-pixel rendered pages inspected; equations, figure, references and margins clear.',
                   initial_issue='3.03pt and 1.03pt paragraph overflows fixed; initial log preserved.'),
        development_repairs=dict(
            initial_exact_session=73242,initial_status='exit 1: all 2520 comparisons agreed, but coverage had two nonzero mixed values and zero negative ones; initial checker preserved.',
            completed_exact_session=1652,completed_status='exit 0; targeted supplemental words added and separately saved.',
            plot_issue='Initial numeric run ended on a matplotlib mathcal label syntax error; braces repaired and full final run passed.'),
        hashes={p.relative_to(ROOT).as_posix():sha(p) for p in files},
        local_render_hashes={p.relative_to(ROOT).as_posix():sha(p) for p in renders},
        build_pass_hashes={f'pass{i}.txt':sha(build/f'pass{i}.txt') for i in [1,2]})
    with target.open('x',encoding='utf-8') as f:json.dump(payload,f,indent=2);f.write('\n')
    print(json.dumps(dict(status='recorded locally',files=len(files),prior_scientific_files=[len(x['unchanged_verified']) for x in prior],pages=8,git_commit=False,git_push=False)))

if __name__=='__main__':main()
