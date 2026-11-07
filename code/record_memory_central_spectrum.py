"""Save central-spectrum provenance and PDF checks; no Git mutation."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,platform,re,subprocess
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-central-spectrum'
BUILD=ROOT/'tmp/pdfs/susskind-memory-central-spectrum-build'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    target=OUT/'local-evidence.json'
    if target.exists():raise FileExistsError(target)
    records={}
    for name in ['trace-ratio','force-bridge','spectrum-assembly']:
        d=json.loads((OUT/(name+'.json')).read_text());assert d['status']=='passed'
        for n,h in d['input_sha256'].items():assert sha(ROOT/n)==h,n
        records[name]=d
    assert len(records['trace-ratio']['records'])==324
    assert len(records['trace-ratio']['finite_exact_checks'])==24
    assert records['force-bridge']['finite_identities']==108
    assert records['spectrum-assembly']['reused_reverse_templates']==128
    parentpath=ROOT/'results/susskind-memory-central-gram/local-evidence.json'
    parent=json.loads(parentpath.read_text());unchanged=[]
    for n,h in parent['hashes'].items():
        if n in ['README.md','PLAN.md','NEXT.md']:continue
        assert sha(ROOT/n)==h,n
        unchanged.append(n)
    log=(OUT/'main.log').read_text(errors='replace')
    assert re.search(r'6\s+pages',log)
    assert not re.search(r'Overfull|undefined|Rerun to get|Fatal error',log)
    pdf=ROOT/'output/pdf/de-sitter-susskind-memory-central-spectrum.pdf'
    assert sha(pdf)==sha(BUILD/'main.pdf')
    info=subprocess.check_output(['pdfinfo',str(pdf)],text=True,errors='replace')
    assert re.search(r'Pages:\s+6\b',info)
    renders=[BUILD/f'page-{i}.png' for i in range(1,7)];assert all(p.exists() for p in renders)
    head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    assert head==parent['current_git_head']
    names=['code/memory_central_trace.py','code/memory_central_force_bridge.py','code/memory_central_spectrum.py',
           'code/build_memory_central_spectrum.ps1','code/record_memory_central_spectrum.py',
           'manuscript/susskind-memory-central-spectrum/main.tex','output/pdf/de-sitter-susskind-memory-central-spectrum.pdf',
           'notes/susskind-memory-central-spectrum.md','literature/susskind-memory-central-spectrum-sources.md',
           'journal/2026-11-06-susskind-memory-central-spectrum.md','README.md','PLAN.md','NEXT.md']
    paths=[ROOT/n for n in names]+sorted(p for p in OUT.rglob('*') if p.is_file())
    result=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),status='Saved locally; no publication',
        scientific_status='Compact central window: quantitative fixed marked-graph trace bridge, general-charge orientation identity, force matrix measure and actual projected W2 core spectrum V2 mu_q2. Raw central U2 descendants/counterterm and F interfaces remain to be extended, followed by charge tails, finite-memory and the bulk dictionary.',
        current_git_head=head,git_commit_created=False,git_push_performed=False,archive_created=False,
        parent_manifest_sha256=sha(parentpath),previous_stage_scientific_files_unchanged=unchanged,
        checks=dict(exact_trace_ratios=24,central_trace_cases=324,precisions=[90,130],
                    CAR_systems=6,CAR_inputs=18,exact_force_identities=108,reused_reverse_templates=128,reused_reverse_words=6072),
        environment=dict(python=platform.python_version(),calculations=[d['environment'] for d in records.values()]),
        build=dict(command='powershell -NoProfile -ExecutionPolicy Bypass -File code/build_memory_central_spectrum.ps1',
                   compiler=log.splitlines()[0],passes_in_final_build=2,final_pages=6,
                   render_command='pdftoppm -r 95 -png output/pdf/de-sitter-susskind-memory-central-spectrum.pdf tmp/pdfs/susskind-memory-central-spectrum-build/page',
                   visual_review='All six final page PNGs viewed at 95 dpi; checked formulas, references and margins.'),
        failures=[],evidence_scope='Analytic balanced-occupation and exponential-overlap bounds transfer existing fixed marked graph results. Exact general-charge CAR supports the strong orientation bridge. Reused endpoint certificates retain their original scope; no new diagonal kernel scan was needed. L2 replacement transfers bounded spectral tests; unbounded moments require uniform integrability.',
        process_state='All calculations, compilation and rendering terminal; no active calculation.',
        hashes={p.relative_to(ROOT).as_posix():sha(p) for p in paths},
        render_hashes={p.relative_to(ROOT).as_posix():sha(p) for p in renders},
        build_pass_hashes={f'pass{i}.txt':sha(BUILD/f'pass{i}.txt') for i in [1,2]})
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status='saved locally',pages=6,files=len(paths),previous_scientific_files_unchanged=len(unchanged),git_head=head)))

if __name__=='__main__':main()
