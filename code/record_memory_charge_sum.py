"""Save local provenance for charge summation; no Git mutation."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,platform,re,subprocess
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-charge-sum'
BUILD=ROOT/'tmp/pdfs/susskind-memory-charge-sum-build'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    target=OUT/'local-evidence.json'
    if target.exists():raise FileExistsError(target)
    records={}
    for name in ['global-sector-checks','growing-window-checks','charge-sum-assembly']:
        data=json.loads((OUT/(name+'.json')).read_text());assert data['status']=='passed'
        for n,h in data['input_sha256'].items():assert sha(ROOT/n)==h,n
        records[name]=data
    assert records['global-sector-checks']['exact_variance_ratio_cases']==6766
    assert len(records['global-sector-checks']['exact_CAR_sectors'])==15
    assert len(records['growing-window-checks']['records'])==12
    parentpath=ROOT/'results/susskind-memory-central-raw/local-evidence.json'
    parent=json.loads(parentpath.read_text());unchanged=[]
    for n,h in parent['hashes'].items():
        if n in ['README.md','PLAN.md','NEXT.md']:continue
        assert sha(ROOT/n)==h,n
        unchanged.append(n)
    log=(OUT/'main.log').read_text(errors='replace')
    assert re.search(r'6\s+pages',log)
    assert not re.search(r'Overfull|undefined|Rerun to get|Fatal error',log)
    pdf=ROOT/'output/pdf/de-sitter-susskind-memory-charge-sum.pdf'
    assert sha(pdf)==sha(BUILD/'main.pdf')
    info=subprocess.check_output(['pdfinfo',str(pdf)],text=True,errors='replace')
    assert re.search(r'Pages:\s+6\b',info)
    renders=[BUILD/f'page-{i}.png' for i in range(1,7)];assert all(p.exists() for p in renders)
    head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    assert head==parent['current_git_head']
    names=['code/memory_charge_global_checks.py','code/memory_charge_growing_window.py',
           'code/assemble_memory_charge_sum.py','code/build_memory_charge_sum.ps1','code/record_memory_charge_sum.py',
           'manuscript/susskind-memory-charge-sum/main.tex','output/pdf/de-sitter-susskind-memory-charge-sum.pdf',
           'notes/susskind-memory-charge-sum.md','literature/susskind-memory-charge-sum-sources.md',
           'journal/2026-11-06-susskind-memory-charge-sum.md','README.md','PLAN.md','NEXT.md']
    links=0
    for n in ['README.md','NEXT.md']:
        text=(ROOT/n).read_text(encoding='utf-8')
        text=text.split('以下为历史继续位置。',1)[0] if n=='NEXT.md' else text.split('## 当前进展：中央原始记忆系数',1)[0]
        for link in re.findall(r'\]\(([^)]+)\)',text):
            if '://' in link or link.startswith('#'):continue
            assert ((ROOT/n).parent/link.split('#')[0]).exists(),link
            links+=1
    paths=[ROOT/n for n in names]+sorted(p for p in OUT.rglob('*') if p.is_file())
    result=dict(status='Saved locally; no publication',created_at_utc=datetime.now(timezone.utc).isoformat(),
        scientific_status='All-sector definition, logarithmic-window remainder audit and charge UI yield the full beta=mu=0 memory-coefficient spectral matrix and Gaussian charge factorization. Exact energy pinching preserves its three Taylor coefficients. Finite-memory covariance remains open.',
        current_git_head=head,git_commit_created=False,git_push_performed=False,archive_created=False,
        parent_manifest_sha256=sha(parentpath),previous_scientific_files_unchanged=unchanged,
        checks=dict(variance_cases=6766,CAR_sectors=15,growing_window_full_rank_groups=12,
                    binomial_tail_sums=4,precisions=[80,120],saved_digits=65,reused_current_certificates=6,local_links=links),
        build=dict(command='powershell -NoProfile -ExecutionPolicy Bypass -File code/build_memory_charge_sum.ps1',
                   passes_in_final_build=2,final_pages=6,compiler=log.splitlines()[0],
                   render_command='pdftoppm -r 95 -png output/pdf/de-sitter-susskind-memory-charge-sum.pdf tmp/pdfs/susskind-memory-charge-sum-build/page',
                   visual_review='All six final pages viewed at 95 dpi; equations, margins and references checked.'),
        environment=dict(python=platform.python_version(),calculations=[d['environment'] for d in records.values()]),
        failures=['First LaTeX build: undefined mathscr; replaced by mathsf.',
                  'Second build: overfull inline system list; split math fragments. Original logs retained.'],
        evidence_scope='Numerical checks support exact finite identities and retain finite-size differences; the growing-window and full-charge result relies on the explicit analytic remainder audit. No fresh whole-memory covariance calculation.',
        process_state='All calculations, compilation and rendering exited successfully; no new active process handle.',
        hashes={p.relative_to(ROOT).as_posix():sha(p) for p in paths},
        render_hashes={p.relative_to(ROOT).as_posix():sha(p) for p in renders},
        build_pass_hashes={f'pass{i}.txt':sha(BUILD/f'pass{i}.txt') for i in [1,2]})
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status='saved locally',pages=6,files=len(paths),prior_scientific_files_unchanged=len(unchanged))))

if __name__=='__main__':main()
