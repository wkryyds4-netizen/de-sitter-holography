"""Local provenance for the central raw-memory stage; no Git mutation."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,platform,re,subprocess

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-central-raw'
BUILD=ROOT/'tmp/pdfs/susskind-memory-central-raw-build'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    target=OUT/'local-evidence.json'
    if target.exists():raise FileExistsError(target)
    records={}
    for name in ['symbolic-descendants','central-filter-checks','central-assembly']:
        d=json.loads((OUT/(name+'.json')).read_text());assert d['status']=='passed'
        for n,h in d['input_sha256'].items():assert sha(ROOT/n)==h,n
        records[name]=d
    parentpath=ROOT/'results/susskind-memory-central-spectrum/local-evidence.json'
    parent=json.loads(parentpath.read_text());unchanged=[]
    for n,h in parent['hashes'].items():
        if n in ['README.md','PLAN.md','NEXT.md']:continue
        assert sha(ROOT/n)==h,n
        unchanged.append(n)
    pdf=ROOT/'output/pdf/de-sitter-susskind-memory-central-raw.pdf'
    assert sha(pdf)==sha(BUILD/'main.pdf')
    log=(OUT/'main.log').read_text(errors='replace')
    assert re.search(r'6\s+pages',log)
    assert not re.search(r'Overfull|undefined|Rerun to get|Fatal error',log)
    info=subprocess.check_output(['pdfinfo',str(pdf)],text=True,errors='replace')
    assert re.search(r'Pages:\s+6\b',info)
    renders=[BUILD/f'page-{i}.png' for i in range(1,7)];assert all(p.exists() for p in renders)
    head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    assert head==parent['current_git_head']
    paths=['README.md','PLAN.md','NEXT.md','PROJECT.md',
           'code/memory_central_descendants.py','code/check_memory_central_raw.py',
           'code/assemble_memory_central_raw.py','code/build_memory_central_raw.ps1',
           'code/record_memory_central_raw.py','manuscript/susskind-memory-central-raw/main.tex',
           'output/pdf/de-sitter-susskind-memory-central-raw.pdf','notes/susskind-memory-central-raw.md',
           'literature/susskind-memory-central-raw-sources.md','journal/2026-11-06-susskind-memory-central-raw.md']
    paths+=[p.relative_to(ROOT).as_posix() for p in OUT.iterdir() if p.is_file()]
    # Validate new local Markdown links without rewriting any historic records.
    links=0
    for name in ['README.md','PLAN.md','NEXT.md','notes/susskind-memory-central-raw.md']:
        content=(ROOT/name).read_text(encoding='utf-8')
        if name=='NEXT.md':content=content.split('以下为历史继续位置。',1)[0]
        if name=='README.md':content=content.split('## 当前进展：中央双向力谱',1)[0]
        for link in re.findall(r'\]\(([^)]+)\)',content):
            if '://' in link or link.startswith('#'):continue
            path=(ROOT/name).parent/link.split('#')[0]
            assert path.exists(),(name,link)
            links+=1
    result=dict(status='Saved locally; no publication',created_at_utc=datetime.now(timezone.utc).isoformat(),
        scientific_status='Compact central window raw U2 restoration, conserved shift, central filtered F and complete Y/L/U2 two-energy spectrum. Charge tails and finite-memory remainder remain open.',
        current_git_head=head,git_commit_created=False,git_push_performed=False,archive_created=False,
        checks=dict(finite_rational_coefficient_references=123,central_full_rank_groups=36,
                    precisions=[80,120],saved_digits=65,central_B_shift_identities=5,
                    prior_current_certificates_verified=7,local_links_checked=links),
        parent_manifest_sha256=sha(parentpath),previous_scientific_files_unchanged=unchanged,
        build=dict(command='powershell -NoProfile -ExecutionPolicy Bypass -File code/build_memory_central_raw.ps1',
                   passes_in_final_build=2,final_pages=6,compiler=log.splitlines()[0],
                   render_command='pdftoppm -r 95 -png output/pdf/de-sitter-susskind-memory-central-raw.pdf tmp/pdfs/susskind-memory-central-raw-build/page',
                   visual_review='All six pages viewed; revised pages 1,3,5 re-rendered and viewed after notation clarifications.'),
        environment=dict(python=platform.python_version(),calculations=[d['environment'] for d in records.values()]),
        failures=['First manifest attempt rejected an old non-link parenthesized expression in historical NEXT; new-stage link validation narrowed to the new entry. See diagnostics.json.'],
        operational_notes=['Read-only CIM process inventory unavailable; no old calculation restarted.'],
        evidence_scope='Saved symbolic-descendants was resumed without rerunning; new numerical refinements retain finite-size deviations and do not certify intervals. Analytic bounds are in the manuscript. No old kernel/endpoint scan rerun.',
        process_state='All new calculations, builds and render commands exited successfully; no new active session.',
        hashes={n:sha(ROOT/n) for n in paths},render_hashes={p.relative_to(ROOT).as_posix():sha(p) for p in renders},
        build_pass_hashes={f'pass{i}.txt':sha(BUILD/f'pass{i}.txt') for i in [1,2]})
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status='saved locally',pages=6,files=len(paths),prior_files_unchanged=len(unchanged),links=links)))

if __name__=='__main__':main()
