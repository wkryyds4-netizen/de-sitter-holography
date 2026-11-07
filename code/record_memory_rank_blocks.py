"""Preserve the rank-block stage and its exact and asymptotic evidence locally."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,re,subprocess,platform
import sympy as s

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-rank-blocks'
BUILD=ROOT/'tmp/pdfs/susskind-memory-rank-blocks-build'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    target=OUT/'local-evidence.json'
    if target.exists():raise FileExistsError(target)
    parent=json.loads((ROOT/'results/susskind-memory-quadratic-contraction/local-evidence.json').read_text())
    previous=[]
    for stage in [r['stage'] for r in parent['previous_stages']]+['susskind-memory-quadratic-contraction']:
        path=ROOT/f'results/{stage}/local-evidence.json';old=json.loads(path.read_text());kept=[]
        for name,digest in old['hashes'].items():
            if name in ['README.md','PLAN.md','NEXT.md']:continue
            assert sha(ROOT/name)==digest,name
            kept.append(name)
        previous.append(dict(stage=stage,manifest_sha256=sha(path),unchanged_verified=kept))
    names=['rank-blocks','rank-reconstruction','product-rank-moments','deficit-limit-checks','product-covariance-checks']
    records={}
    for name in names:
        r=json.loads((OUT/f'{name}.json').read_text());assert r['status']=='passed'
        for path,digest in r['input_sha256'].items():assert sha(ROOT/path)==digest,path
        records[name]=r
    blocks=records['rank-blocks'];rec=records['rank-reconstruction'];mom=records['product-rank-moments']
    assert len(blocks['rows'])==18 and sum(r['gram_rank']==1 for r in blocks['rows'])==11
    assert sum(r['gram_rank']==0 for r in blocks['rows'])==7
    assert blocks['sum_matches_full_and_traceless_previous_gram']
    assert blocks['traceless_global_residual_norm_squared']=='248945355/14903'
    assert len(rec['filter_table'])==11 and sum(r['operator_identities'] for r in rec['finite_reconstructions'])==8
    lam=s.symbols('lam')
    assert all(s.expand(s.sympify(a)-b)==0 for a,b in zip(mom['ds_scaled_delta_casimir_moments'],
                [2*lam,2*lam+4*lam**2,2*lam+12*lam**2+8*lam**3]))
    assert len(records['deficit-limit-checks']['records'])==8
    assert len(records['product-covariance-checks']['rows'])==2
    log=(OUT/'main.log').read_text(errors='replace')
    assert re.search(r'6\s+pages',log)
    assert not re.search(r'Overfull|undefined|Rerun to get|Fatal error',log)
    pdf=ROOT/'output/pdf/de-sitter-susskind-memory-rank-blocks.pdf'
    assert sha(pdf)==sha(BUILD/'main.pdf')
    renders=[BUILD/f'page-{i}.png' for i in range(1,7)];assert all(p.exists() for p in renders)
    head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    assert head==parent['current_git_head']
    files=['code/memory_rank_blocks.py','code/check_memory_rank_reconstruction.py','code/rank_product_moments.py',
           'code/check_rank_deficit_limit.py','code/check_rank_product_covariance.py',
           'code/build_memory_rank_blocks.ps1','code/record_memory_rank_blocks.py',
           'manuscript/susskind-memory-rank-blocks/main.tex','output/pdf/de-sitter-susskind-memory-rank-blocks.pdf',
           'notes/susskind-memory-rank-blocks.md','literature/susskind-memory-rank-blocks-sources.md',
           'journal/2026-11-06-susskind-memory-rank-blocks.md','README.md','PLAN.md','NEXT.md']
    paths=[ROOT/n for n in files]+sorted(p for p in OUT.iterdir() if p.is_file())
    result=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),status='Saved locally; no publication or global freeze',
        scientific_status='Exact finite rank-filter reconstruction; positive ordered-product output-rank law; three general double-scaled deficit moments and tails. General memory filter amplitudes, weighted residual limit and full Z2/finite-memory/bulk dictionary remain open.',
        previous_stages=previous,current_git_head=head,git_commit_created=False,git_push_performed=False,
        archive_created=False,global_project_audit_performed=False,
        checks=dict(rank_blocks=18,nonzero_rank_one_blocks=11,realization_reconstruction_identities=8,
                    general_symbolic_moments=3,finite_racah_transforms=4,independent_product_coordinate_pairs=20000,
                    high_precision_moment_cases=8),
        environment=dict(python=platform.python_version(),calculation_environments=[r['environment'] for r in records.values()],
                         tensor_command_environment='OPENBLAS_NUM_THREADS=1 for rank blocks and CAR product checks'),
        build=dict(command='powershell -NoProfile -ExecutionPolicy Bypass -File code/build_memory_rank_blocks.ps1',
                   compiler=log.splitlines()[0],passes=2,pages=6,
                   visual_review='All six final pages viewed at 95 dpi; no clipping, overflow, overlapping equations, or missing citations.',
                   render_command='pdftoppm -r 95 -png output/pdf/de-sitter-susskind-memory-rank-blocks.pdf tmp/pdfs/susskind-memory-rank-blocks-build/page'),
        retained_issues=['Initial list advanced-index error; code retained.',
                         'Unexpanded symbolic reconstruction comparison; expanded exact differences all zero; code retained.',
                         'Order-three moment requested in rank range R=2; admissible-domain correction retained and larger cases still verify third order.',
                         'Source PDF text available but screenshot and HTML requests returned cache misses; no visual source-inspection claim.'],
        evidence_scope='Rank block equality is full finite quadratic tensor evidence, not matrix samples alone. Output-rank law is for independent isotropic products; same-J symmetrization is bounded, not identified. Three moments do not prove the full Poisson law. Decimal refinements are not interval certificates.',
        process_state='All computation, compilation and rendering terminal. Initial moment29696 exit1; final rank89424, reconstruction60272, moment7602 exit0.',
        hashes={p.relative_to(ROOT).as_posix():sha(p) for p in paths},
        local_render_hashes={p.relative_to(ROOT).as_posix():sha(p) for p in renders},
        build_pass_hashes={f'pass{i}.txt':sha(BUILD/f'pass{i}.txt') for i in [1,2]})
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status='saved locally',files=len(paths),pages=6,
                         prior_scientific_files=sum(len(r['unchanged_verified']) for r in previous),git_head=head)))


if __name__=='__main__':main()
