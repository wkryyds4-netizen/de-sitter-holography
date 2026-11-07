"""Bind the completed local first-order stage without Git mutations."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,re,subprocess
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-first-order'
BUILD=ROOT/'tmp/pdfs/susskind-first-order-build'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    target=OUT/'local-evidence.json'
    if target.exists():raise FileExistsError(target)
    previous=[]
    for stage in ['susskind-singlet-moments','susskind-singlet-double-scale','susskind-singlet-kernel','susskind-singlet-multimark','susskind-singlet-otoc','susskind-rindler-packet']:
        p=ROOT/f'results/{stage}/local-evidence.json';data=json.loads(p.read_text());kept=[]
        for name,h in data['hashes'].items():
            if name in ['README.md','PLAN.md','NEXT.md']:continue
            assert sha(ROOT/name)==h,name
            kept.append(name)
        previous.append(dict(stage=stage,manifest_sha256=sha(p),unchanged_verified=kept))
    for name in ['exact-moments.json','finite-lambda.json','representation-checks.json']:
        d=json.loads((OUT/name).read_text());assert d['status']=='passed'
        for p,h in d['input_sha256'].items():assert sha(ROOT/p)==h,p
    log=(OUT/'main.log').read_text(errors='replace')
    assert not re.search(r'Overfull|undefined|Rerun to get|Fatal error',log)
    assert '(5 pages,' in log
    pdf=ROOT/'output/pdf/de-sitter-susskind-first-order.pdf'
    assert sha(pdf)==sha(BUILD/'main.pdf')
    renders=[BUILD/f'final-page-{i}.png' for i in range(1,6)]
    assert all(p.exists() for p in renders)
    for i in [1,3]:assert sha(BUILD/f'page-{i}.png')==sha(renders[i-1])
    names=['code/singlet_first_order.py','code/check_singlet_first_order_moments.py','code/singlet_general_fourpoint.py','code/check_general_fourpoint.py','code/run_singlet_first_order.py','code/build_singlet_first_order.ps1','code/record_singlet_first_order.py','manuscript/susskind-first-order/main.tex','output/pdf/de-sitter-susskind-first-order.pdf','notes/susskind-first-order.md','literature/susskind-first-order-sources.md','journal/2026-11-06-susskind-first-order.md','README.md','PLAN.md','NEXT.md']
    files=[ROOT/p for p in names]+sorted(p for p in OUT.iterdir() if p.is_file())
    data=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),status='Local stage saved; not published or globally frozen',
        scientific_status='Sourced first-order fixed-time three-channel kernel applied to derivative singlet. Uniform asymptotic remainder and spatial multipoint dictionary remain open.',
        previous_stages=previous,current_git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        git_commit_created=False,git_push_performed=False,archive_created=False,global_project_audit_performed=False,
        checks=dict(exact_taylor_coefficients=22,laurent_cancellations=72,independent_channels=12,finite_lambda_groups=4),
        build=dict(command='powershell -NoProfile -ExecutionPolicy Bypass -File code/build_singlet_first_order.ps1',compiler=log.splitlines()[0],passes=2,pages=5,visual_review='All five pages inspected at 1200px. Final pages 1/3 hashes match inspected images; final pages 2/4/5 re-inspected. No overflow or undefined references.',repaired='Mathtext plotting label; three paragraph overflows; boxed-equation spacing; orphan bibliography.'),
        process_state='All numerical/build/render processes terminal',hashes={p.relative_to(ROOT).as_posix():sha(p) for p in files},
        local_render_hashes={p.relative_to(ROOT).as_posix():sha(p) for p in renders},build_pass_hashes={f'pass{i}.txt':sha(BUILD/f'pass{i}.txt') for i in [1,2]})
    target.write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(dict(status='saved locally',files=len(files),previous_scientific_files=[len(d['unchanged_verified']) for d in previous],pages=5,git_head=data['current_git_head'])))
if __name__=='__main__':main()
