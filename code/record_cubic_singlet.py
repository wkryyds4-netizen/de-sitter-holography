"""Record this local stage without committing or replacing previous evidence."""
from pathlib import Path
from datetime import datetime,timezone
import json,hashlib,re,subprocess
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-cubic-singlet'
BUILD=ROOT/'tmp/pdfs/susskind-cubic-singlet-build'
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def main():
    target=OUT/'local-evidence.json'
    if target.exists():raise FileExistsError(target)
    stages=['susskind-singlet-moments','susskind-singlet-double-scale','susskind-singlet-kernel',
      'susskind-singlet-multimark','susskind-singlet-otoc','susskind-rindler-packet',
      'susskind-first-order','susskind-translations','susskind-antisymmetric','susskind-two-matter']
    previous=[]
    for stage in stages:
        path=ROOT/f'results/{stage}/local-evidence.json';data=json.loads(path.read_text());kept=[]
        for name,h in data['hashes'].items():
            if name in ['README.md','PLAN.md','NEXT.md']:continue
            assert sha(ROOT/name)==h,name
            kept.append(name)
        previous.append(dict(stage=stage,manifest_sha256=sha(path),unchanged_verified=kept))
    for name in ['study.json','spectrum-check.json','identity-check.json']:
        data=json.loads((OUT/name).read_text())
        assert data['status'] in ['passed','computed']
        for path,h in data['input_sha256'].items():assert sha(ROOT/path)==h,path
    study=json.loads((OUT/'study.json').read_text())
    assert study['car']['mark_subset_checks']==3072
    check=json.loads((OUT/'identity-check.json').read_text())
    assert len(check['exact_CAR'])==12 and len(check['symbolic_certificates'])==6
    assert check['two_arc_matches']==9 and check['vanishing_cross_moments']==18
    log=(OUT/'main.log').read_text(errors='replace')
    assert not re.search(r'Overfull|undefined|Rerun to get|Fatal error',log)
    assert '(5 pages,' in log
    assert sha(ROOT/'output/pdf/de-sitter-susskind-cubic-singlet.pdf')==sha(BUILD/'main.pdf')
    renders=[BUILD/f'page-{i}.png' for i in range(1,6)]
    assert all(p.exists() for p in renders)
    files=[ROOT/p for p in [
      'code/cubic_singlet.py','code/check_cubic_spectrum.py','code/check_cubic_identities.py',
      'code/build_cubic_singlet.ps1','code/record_cubic_singlet.py',
      'manuscript/susskind-cubic-singlet/main.tex','output/pdf/de-sitter-susskind-cubic-singlet.pdf',
      'notes/susskind-cubic-singlet.md','literature/susskind-cubic-singlet-sources.md',
      'journal/2026-11-06-susskind-cubic-singlet.md','README.md','PLAN.md','NEXT.md']]
    files+=sorted(p for p in OUT.iterdir() if p.is_file())
    head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    assert head=='714d038bb605cfcffdc351e4dcc62aa95237d25f'
    data=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),
      status='Local stage saved; no publication or global freeze',
      scientific_status='Positive independent cubic Gaussian-chaos singlet and exact moments through four external Hamiltonians. All-order q^3 kernel, raw linear chaos, bulk masses and complete duality remain open.',
      previous_stages=previous,current_git_head=head,git_commit_created=False,git_push_performed=False,
      archive_created=False,global_project_audit_performed=False,
      checks=dict(exact_integer_complex_CAR=12,symbolic_certificates=6,finite_CAR_mark_comparisons=3072,
        two_arc_moment_matches=9,vanishing_cross_moments=18),
      build=dict(command='powershell -NoProfile -ExecutionPolicy Bypass -File code/build_cubic_singlet.ps1',
        compiler=log.splitlines()[0],passes=2,pages=5,
        render_command='pdftoppm -scale-to 1100 -png output/pdf/de-sitter-susskind-cubic-singlet.pdf tmp/pdfs/susskind-cubic-singlet-build/page',
        visual_review='All five final pages visually inspected; no clipping, overlap or unresolved citations.'),
      process_state='All calculation, compilation and rendering sessions terminal',
      hashes={p.relative_to(ROOT).as_posix():sha(p) for p in files},
      local_render_hashes={p.relative_to(ROOT).as_posix():sha(p) for p in renders},
      build_pass_hashes={f'pass{i}.txt':sha(BUILD/f'pass{i}.txt') for i in [1,2]})
    target.write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(dict(status='saved locally',files=len(files),prior_scientific_files=sum(len(p['unchanged_verified']) for p in previous),pages=5,git_head=head)))
if __name__=='__main__':main()
