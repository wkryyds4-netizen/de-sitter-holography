"""Save hashes and scope for the local issue review; no Git writes or external messages."""
from pathlib import Path
import json,hashlib,datetime,re
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/issue-1-review'
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
if __name__=='__main__':
    dest=OUT/'local-evidence.json'
    if dest.exists(): raise FileExistsError(dest)
    result=json.loads((OUT/'independent-assembly.json').read_text())
    for name,h in result['input_hashes'].items(): assert sha(ROOT/name)==h,name
    memory_manifest=ROOT/'results/susskind-memory-return-expansion/local-evidence.json'
    memory=json.loads(memory_manifest.read_text(encoding='utf-8'))
    for name,h in memory['hashes'].items(): assert sha(ROOT/name)==h,name
    log=(OUT/'main.log').read_text(errors='replace')
    assert re.search(r'\(3 pages,',log)
    assert not re.search(r'Overfull|undefined|Rerun to get|Fatal error',log)
    paths=[ROOT/n for n in ['code/review_issue1_evidence.py','code/build_issue1_review.ps1','code/record_issue1_review.py',
        'manuscript/issue-1-review/main.tex','output/pdf/de-sitter-issue-1-review.pdf',
        'notes/issue-1-review.md','literature/issue-1-review-sources.md','journal/2026-11-06-issue-1-review.md',
        'README.md','PLAN.md','NEXT.md']]
    paths+=list((ROOT/'sources/external-reviews/issue-1').glob('*'))
    paths+=[OUT/'independent-assembly.json',OUT/'main.log',OUT/'attachment-text.txt',OUT/'diagnostics/first-layout.log']
    render=ROOT/'tmp/pdfs/issue-1-review-build'
    pngs=[render/f'page-{i}.png' for i in range(1,4)]
    external=[OUT/f'attachment-page-{i}.png' for i in range(1,9)]
    assert all(p.exists() for p in pngs+external)
    evidence=dict(created_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        status='Local initial review complete; overall research goal active',
        external_issue=result['external_issue'],package_hash_checks=43,upstream_hash_checks=18,
        exact_rational_assembly_passed=True,full_cell_reintegration_performed=False,
        recorded_external_replay_failed_checks=9890,external_failure_status_preserved=True,
        no_external_scripts_executed=True,no_external_reply_sent=True,
        git_commit_created=False,git_push_performed=False,archive_created=False,
        scientific_scope='Analytic objections checked; final arithmetic independently reconstructed from supplied integration enclosures. Neither a full independent cubature run nor a completed positive-curvature construction.',
        memory_stage_manifest_sha256=sha(memory_manifest),memory_stage_files_unchanged=len(memory['hashes']),
        build=dict(command='powershell -NoProfile -ExecutionPolicy Bypass -File code/build_issue1_review.ps1',
                   pages=3,passes=2,render_dpi=95,all_pages_visually_inspected=True,
                   external_attachment_pages_read_and_inspected=8,external_attachment_render_dpi=80),
        failures=['Initial network/ZIP timeouts resolved by terminal-state checked retry and resume.',
                  'Initial console decoding failed under gbk; read succeeded with UTF-8 stdout.',
                  'First layout overfull commit hash; preserved first-layout.log; final clean.',
                  'One multi-file patch failed on README context and made no changes; re-applied after inspection.'],
        process_state='All downloads, calculations, builds and renders terminal.',
        hashes={p.relative_to(ROOT).as_posix():sha(p) for p in paths},
        render_hashes={p.relative_to(ROOT).as_posix():sha(p) for p in pngs+external})
    dest.write_text(json.dumps(evidence,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status='passed',pages=3,memory_files_unchanged=len(memory['hashes']),saved=str(dest))))
