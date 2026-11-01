"""Bind the heat cutoff certificate and uncertified matrix study, once."""
from pathlib import Path
import json, hashlib, zipfile
from datetime import datetime, timezone
from audit_heat_matrix import audit

def main():
    root = Path(__file__).resolve().parents[1]
    base = root/'results/stage4/heat-matrix'
    if (base/'evidence.json').exists():
        raise RuntimeError('Already bound')
    sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    priorpath = root/'results/stage4/paired-kernel/evidence.json'
    prior = json.loads(priorpath.read_text())
    for p, h in prior['hashes'].items():
        assert sha(root/p) == h, p
    assert len(prior['hashes']) == 18
    known = {r['name']: r['sha256'] for r in json.loads(
        (root/'results/stage4/paired-kernel/preservation.json').read_text())['old_archives']}
    archives = []
    for p in sorted((root/'snapshots').glob('*.zip')):
        h = sha(p)
        side = p.with_suffix('.sha256')
        assert side.exists() or p.name in known, p
        if side.exists():
            assert side.read_text().split()[0] == h, p
        if p.name in known:
            assert h == known[p.name], p
        archives.append(dict(name=p.name, sha256=h))
    assert len(archives) == 57 and len(known) == 56
    same = 0
    changed = []
    allowed = {'README.md', 'PLAN.md', 'NEXT.md', 'journal/2026-11-01.md',
               'code/project_audit.py', 'results/project-audit.json', 'results/project-files.json'}
    with zipfile.ZipFile(root/'snapshots/2026-11-01-paired-kernel.zip') as z:
        for p in z.namelist():
            if p.endswith('/'):
                continue
            if (root/p).is_file() and sha(root/p) == hashlib.sha256(z.read(p)).hexdigest():
                same += 1
            else:
                assert p in allowed, p
                changed.append(p)
    (base/'preservation.json').write_text(json.dumps(dict(status='passed',
        old_files_unchanged=same, changed_entries=changed,
        prior_bound_files=len(prior['hashes']), old_archives=archives), indent=2))
    sources = {
        'sources/kim-scarinci-2112.13329.pdf': 'https://arxiv.org/pdf/2112.13329v2',
        'sources/garoufalidis-kashaev-1706.08132v2.pdf': 'https://arxiv.org/pdf/1706.08132v2'}
    (base/'source-record.json').write_text(json.dumps(dict(
        sources=[dict(file=p, url=u, sha256=sha(root/p)) for p,u in sources.items()],
        prior_evidence={priorpath.relative_to(root).as_posix(): sha(priorpath)},
        online_references=[
            dict(url='https://dlmf.nist.gov/7.8', reading_scope='Mills ratio and Gaussian tail bounds'),
            dict(url='https://mpmath.readthedocs.io/en/1.3.0/contexts.html',
                 reading_scope='Interval inclusion semantics and experimental status')],
        reading_scope='Cached KS product conventions and flat control; GK Appendix B Fourier framework. Online pages have no content-hash certification. No complete paired-pentagon proof claimed.'), indent=2))
    replay = audit(root)
    assert replay['status'] == 'passed' and replay['checks'] == 188, replay
    paths = ['code/'+p for p in ['heat_matrix.py', 'run_heat_matrix.py', 'direct_heat_left.py',
        'heat_cutoff_bound.py', 'rational_heat_certificate.py', 'audit_heat_matrix.py',
        'build_heat_matrix.ps1', 'record_heat_matrix.py']]
    paths += ['manuscript/heat-matrix/main.tex', 'output/pdf/de-sitter-heat-matrix.pdf',
              'notes/stage4-heat-matrix.md', 'literature/stage4-heat-matrix-sources.md']
    paths += [p.relative_to(root).as_posix() for p in sorted(base.iterdir()) if p.is_file()]
    evidence = dict(created_at_utc=datetime.now(timezone.utc).isoformat(),
        status='diagnostic stage bound; project audit follows', scientific_status=replay['scientific_status'],
        pages=6, matrix_records=37, direct_left_records=12, interval_evaluations=6,
        rational_certificates=3, saved_arithmetic_replay_checks=188,
        compilers=['pdfLaTeX twice exit 0', 'desktop success'],
        visual_review=dict(pages=list(range(1,7)), longest_edge_pixels=1100, status='all inspected'),
        hashes={p:sha(root/p) for p in paths})
    (base/'evidence.json').write_text(json.dumps(evidence, indent=2))
    print(json.dumps(dict(bound_files=len(paths), old_files_unchanged=same,
        old_archives=len(archives), scientific_status=evidence['scientific_status'])))

if __name__ == '__main__':
    main()
