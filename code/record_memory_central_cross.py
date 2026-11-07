"""Record local evidence without Git mutations or publication."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,platform,re,subprocess

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-central-cross'
BUILD=ROOT/'tmp/pdfs/susskind-memory-central-cross-build'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    target=OUT/'local-evidence.json'
    if target.exists():raise FileExistsError(target)
    records={}
    for name in ['car-exploration','charged-casimir','transport','cross-limit']:
        data=json.loads((OUT/(name+'.json')).read_text())
        assert data['status']==('exploratory candidate check' if name=='car-exploration' else 'passed')
        for n,h in data['input_sha256'].items():assert sha(ROOT/n)==h,n
        records[name]=data
    assert all(r['candidate_XA_difference']=='0' for r in records['car-exploration']['rows'])
    assert records['charged-casimir']['exact_identities']==60
    assert len(records['cross-limit']['records'])==36
    assert records['cross-limit']['contiguous_identities']==624
    assert records['transport']['finite']['ordered_gaussian_pairs']==1296
    parentpath=ROOT/'results/susskind-memory-central-norm/local-evidence.json'
    parent=json.loads(parentpath.read_text());unchanged=[]
    for n,h in parent['hashes'].items():
        if n in ['README.md','PLAN.md','NEXT.md']:continue
        assert sha(ROOT/n)==h,n
        unchanged.append(n)
    log=(OUT/'main.log').read_text(errors='replace')
    assert re.search(r'6\s+pages',log)
    assert not re.search(r'Overfull|undefined|Rerun to get|Fatal error',log)
    pdf=ROOT/'output/pdf/de-sitter-susskind-memory-central-cross.pdf'
    assert sha(pdf)==sha(BUILD/'main.pdf')
    info=subprocess.check_output(['pdfinfo',str(pdf)],text=True,errors='replace')
    assert re.search(r'Pages:\s+6\b',info)
    renders=[BUILD/f'page-{i}.png' for i in range(1,7)]
    assert all(p.exists() for p in renders)
    head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    assert head==parent['current_git_head']
    names=['code/explore_memory_central_cross.py','code/memory_central_cross.py',
           'code/check_central_charged_casimir.py','code/check_memory_central_transport.py',
           'code/check_memory_central_cross_limit.py','code/build_memory_central_cross.ps1',
           'code/record_memory_central_cross.py','manuscript/susskind-memory-central-cross/main.tex',
           'output/pdf/de-sitter-susskind-memory-central-cross.pdf',
           'notes/susskind-memory-central-cross.md','literature/susskind-memory-central-cross-sources.md',
           'journal/2026-11-06-susskind-memory-central-cross.md','README.md','PLAN.md','NEXT.md']
    paths=[ROOT/n for n in names]+sorted(p for p in OUT.rglob('*') if p.is_file())
    result=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),status='Saved locally; no publication',
        scientific_status='Even N, compact central window: adjacent product transport O(1/p), output-deficit transport O(p^-2/3), and B/L_input/S Gram proved from prior norm/tail results. General-charge S--O remains conditional on the exchanged gradient recoupling lemma, despite finite checks; no completed four-field Gram or dynamical extension claimed.',
        current_git_head=head,git_commit_created=False,git_push_performed=False,archive_created=False,
        parent_manifest_sha256=sha(parentpath),previous_stage_scientific_files_unchanged=unchanged,
        checks=dict(pure_rank_pairs=15,exact_casimir_identities=60,ordered_CAR_pairs=66531,
                    symmetric_CAR_pairs=31466,transport_coordinate_pairs=1296,
                    contiguous_identities=624,central_all_rank_cases=36,precisions=[80,120],stored_significant_digits=65),
        environment=dict(python=platform.python_version(),calculations=[d['environment'] for d in records.values()]),
        build=dict(command='powershell -NoProfile -ExecutionPolicy Bypass -File code/build_memory_central_cross.ps1',
                   compiler=log.splitlines()[0],passes_in_final_build=2,final_pages=6,
                   render_command='pdftoppm -r 95 -png output/pdf/de-sitter-susskind-memory-central-cross.pdf tmp/pdfs/susskind-memory-central-cross-build/page',
                   visual_review='All six final page PNGs viewed at 95 dpi; no clipping, overlaps or unreadable formula labels.'),
        failures=[dict(stage='initial LaTeX build',reason='Single rather than double backslashes at align row endings produced overflow and multiple labels. First source/log retained in diagnostics; corrected final double build passed.',resolved=True)],
        evidence_scope='Transport proof relies on already established common-rank removal isometry, coherent A weighted norm, factorial input tails and positive sixth output moment. Finite CAR and all-rank sums are supporting checks. Candidate Casimir cross is not upgraded to a theorem.',
        process_state='All calculations, compilation and rendering terminal; no active calculation.',
        hashes={p.relative_to(ROOT).as_posix():sha(p) for p in paths},
        render_hashes={p.relative_to(ROOT).as_posix():sha(p) for p in renders},
        build_pass_hashes={f'pass{i}.txt':sha(BUILD/f'pass{i}.txt') for i in [1,2]})
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status='saved locally',pages=6,files=len(paths),
                         previous_scientific_files_unchanged=len(unchanged),git_head=head)))

if __name__=='__main__':main()
