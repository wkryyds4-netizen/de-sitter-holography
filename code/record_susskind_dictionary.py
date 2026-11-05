"""Record the reviewed dictionary stage without changing prior scientific inputs."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import re
import subprocess
from audit_complex_syk_dictionary import audit
from check_singlet_charge_bridge import verify

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-dictionary'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def save(name,data):
    with (OUT/name).open('x',encoding='utf-8') as f:
        json.dump(data,f,ensure_ascii=False,indent=2);f.write('\n')


def main():
    for n in ['preservation.json','build-review.json','source-record.json','evidence.json']:
        if (OUT/n).exists():
            raise FileExistsError(OUT/n)
    assert audit()['status']=='passed'
    assert verify()['passed']
    old=json.loads((ROOT/'results/project-files.json').read_text(encoding='utf-8'))
    mutable={'README.md','PLAN.md','NEXT.md','PROJECT.md','code/project_audit.py',
             'results/project-audit.json','results/project-files.json'}
    checked=[]
    for row in old:
        name=row['file']
        if name in mutable or name.startswith('journal/'):
            continue
        assert sha(ROOT/name)==row['sha256'],name
        checked.append(name)
    now=datetime.now(timezone.utc).isoformat()
    save('preservation.json',dict(created_at_utc=now,prior_inventory_sha256=sha(ROOT/'results/project-files.json'),
        prior_files_checked=len(checked),all_preserved=True,
        scope='Previous inventory scientific files checked bytewise; existing ZIP archives left untouched, not rehashed in this stage.'))
    log=(OUT/'main.log').read_text(encoding='utf-8',errors='replace')
    assert not re.search(r'undefined|Rerun to get|Overfull|Fatal error',log)
    assert '(6 pages,' in log
    pages=[ROOT/f'tmp/pdfs/susskind-dictionary-build/page-{i}.png' for i in range(1,7)]
    assert all(p.exists() for p in pages)
    version=subprocess.run(['pdflatex','--version'],capture_output=True,check=True).stdout.decode('utf-8',errors='replace').splitlines()[0]
    save('build-review.json',dict(created_at_utc=now,pages=6,latex_passes=2,compiler=version,
        visual_review='All six final 1200-pixel PNG pages inspected; no clipping, overlap or formula overflow.',
        initial_failure='Long script-name paragraph overflow; initial log retained, final text repaired and rebuilt.',
        build_command='powershell -File code/build_susskind_dictionary.ps1',
        local_render_sha256={p.relative_to(ROOT).as_posix():sha(p) for p in pages}))
    ledger=ROOT/'literature/susskind-dictionary-sources.md'
    save('source-record.json',dict(created_at_utc=now,ledger=ledger.relative_to(ROOT).as_posix(),
        urls=sorted(set(re.findall(r'https://arxiv\.org/html/[^)\s]+',ledger.read_text(encoding='utf-8')))),
        reading_scope='Exact inspected passages and unread follow-ups listed in ledger; no new full-text cache.'))
    files=[ROOT/p for p in [
        'code/complex_syk_dictionary.py','code/run_complex_syk_dictionary.py',
        'code/audit_complex_syk_dictionary.py','code/check_singlet_charge_bridge.py',
        'code/build_susskind_dictionary.ps1','code/record_susskind_dictionary.py',
        'manuscript/susskind-dictionary/main.tex','output/pdf/de-sitter-susskind-dictionary.pdf',
        'notes/susskind-dictionary.md','literature/susskind-dictionary-sources.md']]
    files+=sorted(p for p in OUT.iterdir() if p.is_file() and p.name!='evidence.json')
    save('evidence.json',dict(created_at_utc=now,status='recorded',
        scientific_status='Exact finite-system identities, literature normalization map and finite-N calibration; large-N singlet spectrum remains open.',
        hashes={p.relative_to(ROOT).as_posix():sha(p) for p in files}))
    print(json.dumps(dict(bound_files=len(files),old_files_preserved=len(checked),pages=6)))


if __name__=='__main__':
    main()
