"""Record a local, uncommitted research increment without a global freeze."""
from pathlib import Path
from fractions import Fraction
from datetime import datetime, timezone
import hashlib
import json
import re
import subprocess

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-singlet-moments'


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    target=OUT/'local-evidence.json'
    if target.exists():raise FileExistsError(target)
    data=json.loads((OUT/'moments.json').read_text(encoding='utf-8'))
    assert data['status']=='passed'
    for path,digest in data['input_sha256'].items():assert sha(ROOT/path)==digest,path
    mainrow=data['records'][-1]
    assert Fraction(mainrow['string_mu2_Js1'])==Fraction(15,8)
    assert sum(map(Fraction,mainrow['independent_15_pairings']))==120
    assert max(data['independent_cubature']['scaled_residuals'])<1e-12
    log=(OUT/'main.log').read_text(encoding='utf-8',errors='replace')
    assert not re.search(r'Overfull|undefined|Fatal error|Rerun to get',log)
    assert '(5 pages,' in log
    renders=[ROOT/f'tmp/pdfs/susskind-singlet-moments-build/page-{i}.png' for i in range(1,6)]
    assert all(p.exists() for p in renders)
    # Read-only Git checks: no commit, index mutation, fetch or push.
    changed=subprocess.check_output(['git','diff','--name-only'],cwd=ROOT,text=True).splitlines()
    assert set(changed)<={'README.md','PLAN.md','NEXT.md'},changed
    files=[ROOT/p for p in ['code/singlet_wick_moments.py','code/run_singlet_wick_moments.py',
           'code/build_singlet_wick_moments.ps1','code/record_singlet_wick_local.py',
           'manuscript/susskind-singlet-moments/main.tex',
           'output/pdf/de-sitter-susskind-singlet-moments.pdf',
           'notes/susskind-singlet-moments.md','literature/susskind-singlet-moments-sources.md']]
    files+=sorted(p for p in OUT.iterdir() if p.is_file())
    payload=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),
        status='local increment recorded; not frozen or published',
        scientific_status='Finite-N Wick formulas and symmetry reduction; general combinatorial and double-scaled evaluation still pending.',
        current_git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        tracked_edits=changed,git_commit_created=False,git_push_performed=False,
        checks='Saved code hashes, exact moment arithmetic, cubature residuals, final build log and tracked-change scope; no global project audit or archive rehash.',
        build=dict(command='powershell -NoProfile -ExecutionPolicy Bypass -File code/build_singlet_wick_moments.ps1',
            compiler=subprocess.check_output(['pdflatex','--version'],text=True,errors='replace').splitlines()[0],
            passes=2,pages=5,visual_review='All five 1200-pixel PNG pages inspected; no clipping or overlap.',
            initial_issue='Long path overfull paragraph; initial log retained and final layout repaired.'),
        hashes={p.relative_to(ROOT).as_posix():sha(p) for p in files},
        local_render_hashes={p.relative_to(ROOT).as_posix():sha(p) for p in renders})
    with target.open('x',encoding='utf-8') as f:json.dump(payload,f,indent=2);f.write('\n')
    print(json.dumps(dict(status='recorded locally',files=len(files),pages=5,git_commit=False,git_push=False)))


if __name__=='__main__':main()
