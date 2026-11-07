"""Bind OTOC stage, sources and rendered PDF; preserve earlier science."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,re,subprocess

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-singlet-otoc'
BUILD=ROOT/'tmp/pdfs/susskind-singlet-otoc-build'


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    target=OUT/'local-evidence.json'
    if target.exists():raise FileExistsError(target)
    previous=[]
    for stage in ['susskind-singlet-moments','susskind-singlet-double-scale','susskind-singlet-kernel','susskind-singlet-multimark']:
        p=ROOT/f'results/{stage}/local-evidence.json';d=json.loads(p.read_text());kept=[];navigation=[]
        for name,h in d['hashes'].items():
            if name in ['README.md','PLAN.md','NEXT.md']:navigation.append(name)
            else:
                assert sha(ROOT/name)==h,name
                kept.append(name)
        previous.append(dict(stage=stage,manifest_sha256=sha(p),unchanged_verified=kept,intentionally_updated_navigation=navigation))
    checks={}
    for name in ['exact-channels.json','otoc-numerics.json','representation-checks.json','error-replay.json']:
        d=json.loads((OUT/name).read_text());assert d['status']=='passed'
        for p,h in d['input_sha256'].items():assert sha(ROOT/p)==h,p
        checks[name]=d
    assert checks['exact-channels.json']['rational_channel_comparisons']==780
    num=checks['otoc-numerics.json']
    assert sum(len(g['records']) for g in num['groups'])==155
    assert sum(len(g['records']) for g in num['scrambling_scaling'])==9
    assert checks['error-replay.json']['arithmetic_records']==164
    assert len(checks['error-replay.json']['bounds'])==8
    log=(OUT/'main.log').read_text(errors='replace')
    assert not re.search(r'Overfull|undefined|Fatal error|Rerun to get',log)
    assert '(6 pages,' in log
    pdf=ROOT/'output/pdf/de-sitter-susskind-singlet-otoc.pdf'
    assert sha(pdf)==sha(BUILD/'main.pdf')
    renders=[BUILD/f'page-{i}.png' for i in range(1,7)]
    assert all(p.exists() for p in renders)
    for i in [1,2]:
        raw=(BUILD/f'pass{i}.txt').read_bytes()
        enc='utf-16' if raw.startswith((b'\xff\xfe',b'\xfe\xff')) else 'utf-8'
        assert '(6 pages,' in raw.decode(enc,errors='replace')
    source=ROOT/'sources/mukhametzhanov-2303.03474v3.pdf'
    assert source.read_bytes().startswith(b'%PDF')
    source_render=BUILD/'mukh-14.png';assert source_render.exists()
    names=[
        'code/singlet_otoc.py','code/check_singlet_otoc.py','code/run_singlet_otoc.py',
        'code/check_otoc_representation.py','code/check_otoc_saved_evidence.py',
        'code/build_singlet_otoc.ps1','code/record_singlet_otoc_local.py',
        'manuscript/susskind-singlet-otoc/main.tex','output/pdf/de-sitter-susskind-singlet-otoc.pdf',
        'notes/susskind-singlet-otoc.md','literature/susskind-singlet-otoc-sources.md',
        'journal/2026-11-06-susskind-singlet-otoc.md','README.md','PLAN.md','NEXT.md']
    files=[ROOT/p for p in names]+sorted(p for p in OUT.iterdir() if p.is_file())
    payload=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),
        status='local research stage recorded, not published or globally frozen',
        scientific_status='Established BINT/BNR kernels applied to the conserved-subtracted singlet, all same-species channels retained. Finite-positive-lambda alternating-time OTOC computed with analytic cutoff bounds. Known scrambling asymptotic mapped into string units. Bulk wavepacket/time dictionary and complete de Sitter duality remain open.',
        previous_stages=previous,current_git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        git_commit_created=False,git_push_performed=False,archive_created=False,global_project_audit_performed=False,
        build=dict(command='powershell -NoProfile -ExecutionPolicy Bypass -File code/build_singlet_otoc.ps1',
            compiler=log.splitlines()[0],passes=2,pages=6,
            visual_review='All six 1200px pages inspected; final amendment to Appendix A re-rendered and inspected. Formulas, plot, tables, footnote and references clear.',
            initial_issue='0.95389pt abstract overflow repaired; initial log saved.'),
        source_cache=dict(path=source.relative_to(ROOT).as_posix(),sha256=sha(source),url='https://arxiv.org/pdf/2303.03474v3',
            inspected_page=14,printed_page=12,render_sha256=sha(source_render),
            acquisition='Initial sandbox socket access denied; authorized public-paper download via elevated tool succeeded. No approval rejection.',
            text_extraction='pdftotext reported an xref reconstruction warning; source equation page verified visually.'),
        process_state='All numerical, checking and build commands terminal; no live process to poll.',
        hashes={p.relative_to(ROOT).as_posix():sha(p) for p in files},
        local_render_hashes={p.relative_to(ROOT).as_posix():sha(p) for p in renders},
        build_pass_hashes={f'pass{i}.txt':sha(BUILD/f'pass{i}.txt') for i in [1,2]})
    target.write_text(json.dumps(payload,indent=2)+'\n')
    print(json.dumps(dict(status='recorded locally',files=len(files),prior_scientific_files=[len(p['unchanged_verified']) for p in previous],
        pages=6,git_head=payload['current_git_head'],git_commit=False,git_push=False)))

if __name__=='__main__':main()
