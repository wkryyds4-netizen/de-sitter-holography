"""Record local research evidence without committing or publishing."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,re,subprocess
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-rindler-packet'
BUILD=ROOT/'tmp/pdfs/susskind-rindler-packet-build'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    target=OUT/'local-evidence.json'
    if target.exists():raise FileExistsError(target)
    previous=[]
    for stage in ['susskind-singlet-moments','susskind-singlet-double-scale','susskind-singlet-kernel','susskind-singlet-multimark','susskind-singlet-otoc']:
        p=ROOT/f'results/{stage}/local-evidence.json';d=json.loads(p.read_text());kept=[];nav=[]
        for name,h in d['hashes'].items():
            if name in ['README.md','PLAN.md','NEXT.md']:nav.append(name)
            else:
                assert sha(ROOT/name)==h,name
                kept.append(name)
        previous.append(dict(stage=stage,manifest_sha256=sha(p),unchanged_verified=kept,intentionally_updated_navigation=nav))
    for filename in ['study.json','independent-checks.json']:
        d=json.loads((OUT/filename).read_text());assert d['status']=='passed'
        for p,h in d['input_sha256'].items():assert sha(ROOT/p)==h,p
    study=json.loads((OUT/'study.json').read_text())
    assert len(study['checks'])==91 and all(c['passed'] for c in study['checks'])
    independent=json.loads((OUT/'independent-checks.json').read_text())
    assert sum(len(independent[k]) for k in ['radial','oscillators','finite_box_norms','kms'])==19
    log=(OUT/'main.log').read_text(errors='replace')
    assert not re.search(r'Overfull|undefined|Rerun to get|Fatal error',log)
    assert '(5 pages,' in log
    pdf=ROOT/'output/pdf/de-sitter-susskind-rindler-packet.pdf'
    assert sha(pdf)==sha(BUILD/'main.pdf')
    renders=[BUILD/f'final-page-{i}.png' for i in range(1,6)]
    assert all(p.exists() for p in renders)
    for i in range(1,5):assert sha(BUILD/f'page-{i}.png')==sha(renders[i-1]),i
    names=['code/rindler_singlet_packet.py','code/run_rindler_singlet_packet.py','code/check_rindler_singlet_packet.py',
        'code/build_rindler_singlet_packet.ps1','code/record_rindler_singlet_packet.py',
        'manuscript/susskind-rindler-packet/main.tex','output/pdf/de-sitter-susskind-rindler-packet.pdf',
        'notes/susskind-rindler-packet.md','literature/susskind-rindler-packet-sources.md',
        'journal/2026-11-06-susskind-rindler-packet.md','README.md','PLAN.md','NEXT.md']
    files=[ROOT/p for p in names]+sorted(p for p in OUT.iterdir() if p.is_file())
    data=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),
        status='Local stage saved; not published or globally frozen',
        scientific_status='Explicit regulated free Rindler midpoint construction via established KL transform; mass/radial dilation freedom and a specific equal-gap fourpoint mismatch. Full multipoint dictionary and de Sitter duality remain open.',
        previous_stages=previous,current_git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        git_commit_created=False,git_push_performed=False,archive_created=False,global_project_audit_performed=False,
        checks=dict(main=91,independent=19),
        build=dict(command='powershell -NoProfile -ExecutionPolicy Bypass -File code/build_rindler_singlet_packet.ps1',
            compiler=log.splitlines()[0],passes=2,pages=5,
            visual_review='All five 1200px pages inspected. Final pages 1-4 pixel hashes identical to inspected renders; final page 5 re-inspected after reference layout adjustment.',
            repaired_issues='Missing sech declaration; path overflow; orphan final bibliography entry. Initial failed compilation log retained.'),
        sources='Primary sources read online; URLs and equation locations in literature note, no new source-PDF cache.',
        process_state='All numerical, build and rendering commands terminal; no running session.',
        hashes={p.relative_to(ROOT).as_posix():sha(p) for p in files},
        local_render_hashes={p.relative_to(ROOT).as_posix():sha(p) for p in renders},
        build_pass_hashes={f'pass{i}.txt':sha(BUILD/f'pass{i}.txt') for i in [1,2]})
    target.write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(dict(status='saved locally',files=len(files),previous_scientific_files=[len(p['unchanged_verified']) for p in previous],
        pages=5,git_head=data['current_git_head'],git_commit=False,git_push=False)))


if __name__=='__main__':main()
