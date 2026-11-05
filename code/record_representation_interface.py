"""Bind the analytic interface note; preserve earlier science and archives."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import platform
import zipfile
from audit_representation_interface import audit


def main():
    root = Path(__file__).resolve().parents[1]
    base = root / 'results/stage4/representation-interface'
    for name in ('evidence.json', 'preservation.json', 'input-audit.json', 'source-record.json'):
        if (base / name).exists():
            raise RuntimeError('Already recorded: ' + name)
    sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    result = audit(root)
    if result['status'] != 'passed':
        raise RuntimeError(result['errors'])
    old = root / 'results/stage4/residual-certificate'
    known = {r['name']: r['sha256'] for r in
             json.loads((old / 'preservation.json').read_text())['old_archives']}
    archives = []
    for path in sorted((root / 'snapshots').glob('*.zip')):
        digest = sha(path)
        side = path.with_suffix('.sha256')
        if side.exists():
            if side.read_text().split()[0] != digest:
                raise RuntimeError('Archive checksum mismatch: ' + path.name)
        elif path.name not in known:
            raise RuntimeError('No prior checksum: ' + path.name)
        if path.name in known and known[path.name] != digest:
            raise RuntimeError('Prior archive changed: ' + path.name)
        archives.append(dict(name=path.name, sha256=digest))
    if len(archives) != 59 or len(known) != 58:
        raise RuntimeError('Unexpected archive count')
    allowed = {'README.md', 'PLAN.md', 'NEXT.md', 'journal/2026-11-01.md',
               'code/project_audit.py', 'results/project-audit.json', 'results/project-files.json'}
    same, changed = 0, []
    with zipfile.ZipFile(root / 'snapshots/2026-11-01-residual-certificate.zip') as archive:
        for name in archive.namelist():
            if name.endswith('/'):
                continue
            path = root / name
            if path.is_file() and sha(path) == hashlib.sha256(archive.read(name)).hexdigest():
                same += 1
            elif name in allowed:
                changed.append(name)
            else:
                raise RuntimeError('Frozen scientific file changed: ' + name)
    def write(name, value):
        with (base / name).open('x', encoding='utf-8') as stream:
            json.dump(value, stream, indent=2, ensure_ascii=False)
            stream.write('\n')
    write('input-audit.json', result)
    write('preservation.json', dict(status='passed', old_files_unchanged=same,
          changed_entries=changed, old_archives=archives))
    sources = [
        dict(url='https://arxiv.org/abs/2112.13329v2', location='Section 6.3 and Proposition 6.15',
             file='sources/kim-scarinci-2112.13329.pdf',
             sha256=sha(root / 'sources/kim-scarinci-2112.13329.pdf')),
        dict(url='https://arxiv.org/html/2112.13962v2', location='Propositions 7.3 and 7.4; proof dependency text'),
        dict(url='https://arxiv.org/pdf/math-ph/0511030v2', location='Sections 3.1-3.2, Theorem 4, printed page 12',
             access='online extracted PDF text; screenshot request failed; no local download'),
        dict(url='https://link.springer.com/article/10.1007/s00220-024-05012-8', location='article landing page; correction search'),
        dict(url='https://sites.google.com/site/hyunkyukimmath/papers', location='paper entries; correction search')]
    write('source-record.json', dict(checked_date='2026-11-05', sources=sources,
          scope='Local passage reading and targeted correction search, not full verification or exhaustive search.'))
    review = json.loads((base / 'build-review.json').read_text())
    if review['pages'] != 4 or review['visual_review']['status'] != 'all inspected':
        raise RuntimeError('Build review incomplete')
    files = ['code/audit_representation_interface.py', 'code/record_representation_interface.py',
             'code/build_representation_interface.ps1', 'manuscript/representation-interface/main.tex',
             'output/pdf/de-sitter-representation-interface.pdf',
             'notes/stage4-representation-interface.md', 'literature/stage4-representation-interface-sources.md']
    files += [p.relative_to(root).as_posix() for p in sorted(base.iterdir()) if p.is_file()]
    write('evidence.json', dict(created_at_utc=datetime.now(timezone.utc).isoformat(),
          status='stage bound', scientific_status='Analytic consequences of the prior computer-assisted bound under unchanged spectral assumptions',
          scope='No independent recertification, no numerical interval radius, no finite-G claim.',
          python_version=platform.python_version(), pages=4,
          hashes={name: sha(root / name) for name in files}))
    print(json.dumps(dict(bound_files=len(files), old_files_unchanged=same,
                          old_archives=len(archives), input_audit=result['status'])))


if __name__ == '__main__':
    main()
