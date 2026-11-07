"""Record local mixed-singlet evidence and verify prior scientific file hashes."""
from pathlib import Path
from datetime import datetime,timezone
import json,hashlib,platform,re,subprocess
import numpy as np
import sympy as s
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-mixed-cubic'
BUILD=ROOT/'tmp/pdfs/susskind-mixed-cubic-build'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    target=OUT/'local-evidence.json';summary_path=OUT/'verification-summary.json'
    if target.exists() or summary_path.exists():raise FileExistsError('Evidence already exists')
    prior=ROOT/'results/susskind-filtered-spectrum/local-evidence.json'
    parent=json.loads(prior.read_text());stages=[r['stage'] for r in parent['previous_stages']]+['susskind-filtered-spectrum']
    previous=[]
    for stage in stages:
        path=ROOT/f'results/{stage}/local-evidence.json';record=json.loads(path.read_text());kept=[]
        for name,h in record['hashes'].items():
            if name in ['README.md','PLAN.md','NEXT.md']:continue
            assert sha(ROOT/name)==h,name;kept.append(name)
        previous.append(dict(stage=stage,manifest_sha256=sha(path),unchanged_verified=kept))
    data={}
    for name in ['moments','cluster-certificate','finite-car','otoc']:
        r=json.loads((OUT/f'{name}.json').read_text());data[name]=r
        assert r['status']==('computed' if name=='moments' else 'passed')
        for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h,p
    moments=data['moments']['records'];clusters=data['cluster-certificate'];car=data['finite-car'];otoc=data['otoc']
    assert len(moments)==12 and all(r['matches_independent_probe'] and r['independent_probe_difference']=='0' for r in moments)
    assert sum(r['pairings'] for r in moments)==15392 and sum(r['minor_checks'] for r in moments)==433964
    assert clusters['templates']==14 and clusters['pairings']==148
    assert car['comparisons']==6656 and car['words']==64 and car['nonzero']==1031 and car['negative']==2
    q,x=s.symbols('q x');limits=[]
    for r in otoc['coefficients']:
        n=r['order'];c=s.sympify(r['coefficient_in_sigma_t'])
        expr=(2*(1-c) if n==0 else -2*c)/x**(n//2+1)
        val=s.limit(expr.subs(q,s.exp(-4*x)),x,0,dir='+')
        limits.append(dict(time_order=n,small_lambda_coefficient=str(val)))
    assert [r['small_lambda_coefficient'] for r in limits]==['48','384','512','4096/15']
    points=[]
    for r in otoc['records']:
        bound=max(b['otoc_truncation_bound'] for b in r['bounds'])
        for p in r['points']:
            points.append(dict(lambda_value=r['lambda_value'],time=p['time'],
                relative_cutoff_change=p['cutoff_difference']/abs(p['otoc']),
                relative_species_swap_change=p['swap_difference']/abs(p['otoc']),
                conservative_relative_truncation_bound=bound/(abs(p['otoc'])-bound)))
    summary=dict(status='passed',exact_mixed_words=len(moments),pairings=15392,minor_checks=433964,
        cluster_templates=14,cluster_pairings=148,finite_car_words=64,finite_car_comparisons=6656,
        fixed_time_taylor_small_lambda=limits,numerical_comparisons=points,
        maximum_absolute_changes={k:max(p.get(k,0) for r in otoc['records'] for p in r['points']) for k in ['cutoff_difference','swap_difference','series_difference']},
        maximum_relative_changes={k:max(p[k] for p in points) for k in ['relative_cutoff_change','relative_species_swap_change','conservative_relative_truncation_bound']},
        error_scope='Truncation analytically bounded in exact arithmetic. Floating point empirical. Positive OTOC sum; signed Taylor cancellations performed exactly before evaluation.',
        input_sha256={f'results/susskind-mixed-cubic/{n}.json':sha(OUT/f'{n}.json') for n in data})
    log=(OUT/'main.log').read_text(errors='replace');assert not re.search(r'Overfull|undefined|Rerun to get|Fatal error',log) and '(6 pages,' in log
    pdf=ROOT/'output/pdf/de-sitter-susskind-mixed-cubic.pdf';assert sha(pdf)==sha(BUILD/'main.pdf')
    renders=[BUILD/f'final-page-{i}.png' for i in range(1,7)];assert all(p.exists() for p in renders)
    for i in [2,3,4,6]:assert sha(BUILD/f'page-{i}.png')==sha(BUILD/f'final-page-{i}.png'),i
    head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip();assert head=='714d038bb605cfcffdc351e4dcc62aa95237d25f'
    summary_path.write_text(json.dumps(summary,indent=2)+'\n')
    names=['code/mixed_cubic_vertices.py','code/mixed_cubic_clusters.py','code/check_mixed_cubic_car.py','code/mixed_cubic_otoc.py',
      'code/build_mixed_cubic.ps1','code/record_mixed_cubic.py','manuscript/susskind-mixed-cubic/main.tex',
      'output/pdf/de-sitter-susskind-mixed-cubic.pdf','notes/susskind-mixed-cubic.md','literature/susskind-mixed-cubic-sources.md',
      'journal/2026-11-06-susskind-mixed-cubic.md','README.md','PLAN.md','NEXT.md']
    files=[ROOT/n for n in names]+sorted(p for p in OUT.iterdir() if p.is_file())
    evidence=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),status='Local stage saved; no publication or global freeze',
      scientific_status='Fixed-degree H/X/Rhat joint multiprobe law from local component factorization and exact templates. Unequal-probe mixed OTOC evaluated. Bulk dictionary and complete de Sitter duality remain open.',
      previous_stages=previous,current_git_head=head,git_commit_created=False,git_push_performed=False,archive_created=False,
      global_project_audit_performed=False,environment=dict(python=platform.python_version(),numpy=np.__version__,sympy=s.__version__),
      build=dict(command='powershell -NoProfile -ExecutionPolicy Bypass -File code/build_mixed_cubic.ps1',compiler=log.splitlines()[0],passes=2,pages=6,
        rendering='pdftoppm -r 95 -png output/pdf/de-sitter-susskind-mixed-cubic.pdf tmp/pdfs/susskind-mixed-cubic-build/final-page',
        visual_review='All six initial pages inspected; modified pages 1 and 5 reinspected; pages 2,3,4,6 pixel-identical'),
      process_state='All computation, compilation and rendering processes terminal',
      procedural_corrections=['One read-only lookup used a nonexistent old script filename','One atomic document patch did not match NEXT heading; no files applied; corrected patch applied successfully','Crossing product notation clarified to unordered chord pairs before final build'],
      hashes={p.relative_to(ROOT).as_posix():sha(p) for p in files},
      local_render_hashes={p.relative_to(ROOT).as_posix():sha(p) for p in renders},
      build_pass_hashes={f'pass{i}.txt':sha(BUILD/f'pass{i}.txt') for i in [1,2]})
    target.write_text(json.dumps(evidence,indent=2)+'\n')
    print(json.dumps(dict(status='saved locally',files=len(files),prior_scientific_files=sum(len(p['unchanged_verified']) for p in previous),
        small_lambda=limits,relative_checks=summary['maximum_relative_changes'],pages=6,git_head=head)))
if __name__=='__main__':main()
