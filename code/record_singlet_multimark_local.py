"""Bind local multi-singlet evidence and verify earlier scientific files."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,re,subprocess,math

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-singlet-multimark'
BUILD=ROOT/'tmp/pdfs/susskind-singlet-multimark-build'


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    target=OUT/'local-evidence.json'
    if target.exists():raise FileExistsError(target)
    prior=[]
    for stage in ['susskind-singlet-moments','susskind-singlet-double-scale','susskind-singlet-kernel']:
        record=ROOT/f'results/{stage}/local-evidence.json';data=json.loads(record.read_text())
        unchanged=[];navigation=[]
        for p,digest in data['hashes'].items():
            if p in ['README.md','PLAN.md','NEXT.md']:navigation.append(p)
            else:
                assert sha(ROOT/p)==digest,p
                unchanged.append(p)
        prior.append(dict(stage=stage,manifest_sha256=sha(record),unchanged_verified=unchanged,intentionally_updated_navigation=navigation))
    results={}
    for name in ['exact-checks.json','finite-car.json','fourpoint-numerics.json']:
        d=json.loads((OUT/name).read_text());assert d['status']=='passed'
        for p,h in d['input_sha256'].items():assert sha(ROOT/p)==h,p
        results[name]=d
    assert results['exact-checks.json']['words']==95
    assert results['exact-checks.json']['local_graphs']==41210
    car=results['finite-car.json']
    assert car['comparisons']==3584 and car['full_mark_nonzero']==42 and car['full_mark_negative']==0
    numeric=results['fourpoint-numerics.json'];replay=0
    assert len(numeric['exact_collapsed_moments'])==25
    for r in numeric['rows']:
        lam=r['lambda_value'];q=math.exp(-4*lam)
        assert abs(r['kappa'][0]-math.expm1(-16*lam))<1e-12
        assert r['kappa'][20]>0
        for g2,g4,k in zip(r['g2'],r['g4'],r['kappa']):
            assert abs(k-((1+q**4)*g4-2*g2*g2))<1e-14
            replay+=1
        assert r['analytic_connected_truncation_bound']<6e-13
    log=(OUT/'main.log').read_text(errors='replace')
    assert not re.search(r'Overfull|undefined|Fatal error|Rerun to get',log)
    assert '(7 pages,' in log
    pdf=ROOT/'output/pdf/de-sitter-susskind-singlet-multimark.pdf'
    assert sha(pdf)==sha(BUILD/'main.pdf')
    renders=[BUILD/f'page-{i}.png' for i in range(1,8)]
    assert all(p.exists() for p in renders)
    for i in [1,2]:
        raw=(BUILD/f'pass{i}.txt').read_bytes()
        encoding='utf-16' if raw.startswith((b'\xff\xfe',b'\xfe\xff')) else 'utf-8'
        assert '(7 pages,' in raw.decode(encoding,errors='replace')
    inputs=[
        'code/singlet_multimark.py','code/check_singlet_multimark.py','code/check_multimark_car.py',
        'code/run_singlet_fourpoint.py','code/build_singlet_multimark.ps1','code/record_singlet_multimark_local.py',
        'manuscript/susskind-singlet-multimark/main.tex','output/pdf/de-sitter-susskind-singlet-multimark.pdf',
        'notes/susskind-singlet-multimark.md','literature/susskind-singlet-multimark-sources.md',
        'journal/2026-11-06-susskind-singlet-multimark.md','README.md','PLAN.md','NEXT.md']
    files=[ROOT/p for p in inputs]+sorted(p for p in OUT.iterdir() if p.is_file())
    payload=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),
        status='local research stage recorded; not published or globally frozen',
        scientific_status='All fixed joint polynomial moments and fixed real-time correlators reduce to mixed-q chords after conserved subtraction. Collapsed connected four-point slice computed. General OTOC evaluation and bulk dictionary remain open.',
        previous_stages=prior,current_git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        git_commit_created=False,git_push_performed=False,archive_created=False,global_project_audit_performed=False,
        numeric_arithmetic_replayed=replay,
        build=dict(command='powershell -NoProfile -ExecutionPolicy Bypass -File code/build_singlet_multimark.ps1',
            compiler=log.splitlines()[0],passes=2,pages=7,
            visual_review='All seven final 1200-pixel pages inspected, including formulas, plot, table, footnote and references; no clipping or overlap.',
            initial_issue='10.50531pt inline-formula paragraph overflow fixed; initial log retained.'),
        development=dict(initial_car_session=51714,initial_car_status='Stopped with exit 1 at factorial-time 12-slot Laplace determinant; no identity failure. Initial source retained.',
            repair='Exact Bareiss determinant, independently checked on 140 integer matrices against Laplace; final 160 words and 3584 CAR comparisons passed.',
            final_process_status='All computational and build processes completed; no pending session.'),
        hashes={p.relative_to(ROOT).as_posix():sha(p) for p in files},
        local_render_hashes={p.relative_to(ROOT).as_posix():sha(p) for p in renders},
        build_pass_hashes={f'pass{i}.txt':sha(BUILD/f'pass{i}.txt') for i in [1,2]})
    target.write_text(json.dumps(payload,indent=2)+'\n')
    print(json.dumps(dict(status='recorded locally',files=len(files),prior_scientific_files=[len(p['unchanged_verified']) for p in prior],pages=7,replayed=replay,git_head=payload['current_git_head'],git_commit=False,git_push=False)))

if __name__=='__main__':main()
