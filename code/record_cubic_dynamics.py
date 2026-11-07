"""Bind this local stage to its calculations, PDF and preceding science."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,re,subprocess
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-cubic-dynamics'
BUILD=ROOT/'tmp/pdfs/susskind-cubic-dynamics-build'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    path=OUT/'local-evidence.json'
    if path.exists():raise FileExistsError(path)
    old=json.loads((ROOT/'results/susskind-cubic-linear/local-evidence.json').read_text())
    stages=[row['stage'] for row in old['previous_stages']]+['susskind-cubic-linear']
    previous=[]
    for stage in stages:
        source=ROOT/f'results/{stage}/local-evidence.json';record=json.loads(source.read_text());kept=[]
        for name,h in record['hashes'].items():
            if name in ['README.md','PLAN.md','NEXT.md']:continue
            assert sha(ROOT/name)==h,name
            kept.append(name)
        previous.append(dict(stage=stage,manifest_sha256=sha(source),unchanged_verified=kept))
    for name in ['validation.json','two-energy-validation.json','double-scale.json']:
        r=json.loads((OUT/name).read_text());assert r['status']=='passed'
        for n,h in r['input_sha256'].items():assert sha(ROOT/n)==h,n
    v=json.loads((OUT/'validation.json').read_text())
    assert v['exact_covariance_inversions']==1515 and v['exact_factorial_moments']==928
    assert len(v['car_records'])==5 and v['car_records'][-1]['direct_M2']=='3240'
    e=json.loads((OUT/'two-energy-validation.json').read_text())
    assert len(e['all_pairings'])==3 and len(e['aggregate_vs_exact'])==6
    assert e['all_pairings'][-1]['B20_pairings']==['3240','180','0']
    d=json.loads((OUT/'double-scale.json').read_text());assert len(d['records'])==12
    log=(OUT/'main.log').read_text(errors='replace')
    assert not re.search(r'Overfull|undefined|Rerun to get|Fatal error',log) and '(6 pages,' in log
    pdf=ROOT/'output/pdf/de-sitter-susskind-cubic-dynamics.pdf'
    assert sha(pdf)==sha(BUILD/'main.pdf')
    renders=[BUILD/f'checked-page-{i}.png' for i in range(1,7)]
    assert all(p.exists() for p in renders)
    for n in [1,2,3,5,6]:assert sha(BUILD/f'final-page-{n}.png')==sha(BUILD/f'checked-page-{n}.png')
    names=['code/cubic_linear_dynamics.py','code/check_cubic_dynamics.py',
      'code/check_cubic_dynamics_moments.py','code/build_cubic_dynamics.ps1','code/record_cubic_dynamics.py',
      'manuscript/susskind-cubic-dynamics/main.tex','output/pdf/de-sitter-susskind-cubic-dynamics.pdf',
      'notes/susskind-cubic-dynamics.md','literature/susskind-cubic-dynamics-sources.md',
      'journal/2026-11-06-susskind-cubic-dynamics.md','README.md','PLAN.md','NEXT.md']
    files=[ROOT/n for n in names]+sorted(p for p in OUT.iterdir() if p.is_file())
    head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    assert head=='714d038bb605cfcffdc351e4dcc62aa95237d25f'
    record=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),status='Local stage saved; no publication or global freeze',
      scientific_status='Exact rank crossing and same-disorder second dynamic moment; DS limit and a uniform fourth-moment bound establish surviving dynamic weight. Complete F kernel, mixed correlators and the bulk duality remain open.',
      previous_stages=previous,current_git_head=head,git_commit_created=False,git_push_performed=False,
      archive_created=False,global_project_audit_performed=False,
      checks=dict(exact_covariance=1515,factorial_moments=928,direct_dynamic_sectors=5,
        direct_two_energy_sectors=3,exact_aggregate_comparisons=6,double_scale_cases=12,precisions=[80,120],
        precision_comparison_scope='65-digit stored normalized moment values; not interval certification'),
      build=dict(command='powershell -NoProfile -ExecutionPolicy Bypass -File code/build_cubic_dynamics.ps1',
        compiler=log.splitlines()[0],passes=2,pages=6,
        visual_review='All six pages inspected; final page 4 re-inspected after arrow notation repair, other five pages pixel-identical.'),
      process_state='All calculation, compilation and rendering processes terminal',
      hashes={p.relative_to(ROOT).as_posix():sha(p) for p in files},
      local_render_hashes={p.relative_to(ROOT).as_posix():sha(p) for p in renders},
      build_pass_hashes={f'pass{i}.txt':sha(BUILD/f'pass{i}.txt') for i in [1,2]})
    path.write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(dict(status='saved locally',files=len(files),prior_scientific_files=sum(len(x['unchanged_verified']) for x in previous),pages=6,git_head=head)))
if __name__=='__main__':main()
