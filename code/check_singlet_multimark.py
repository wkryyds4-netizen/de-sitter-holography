"""Exact tests of the proposed all-polynomial-moment singlet reduction."""
from pathlib import Path
from itertools import combinations
import hashlib,json,platform,time
from singlet_multimark import raw_moment,predicted,probe_moment

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-singlet-multimark'


def words(m,n):
    for positions in combinations(range(m+n),m):
        yield ''.join('X' if i in positions else 'H' for i in range(m+n))


def main():
    OUT.mkdir(exist_ok=True);target=OUT/'exact-checks.json'
    if target.exists():raise FileExistsError(target)
    rows=[];start=time.monotonic()
    for m,n in [(1,0),(1,2),(1,4),(2,0),(2,2),(2,4),(3,0),(3,2),(3,4),(4,0),(4,2),(5,0),(6,0)]:
        for word in words(m,n):
            actual,count=raw_moment(word);expected=predicted(word)
            assert actual==expected,(word,actual,expected)
            rows.append(dict(word=word,graphs=count,probe_q_polynomial={str(k[1]):v for k,v in probe_moment(word).items()},raw_nonzero_terms=len(actual)))
        print(m,n,'passed',round(time.monotonic()-start,2),flush=True)
    inputs=['code/singlet_multimark.py','code/check_singlet_multimark.py','code/singlet_coefficient.py']
    payload=dict(status='passed',scope='exact polynomial identities from local two-site minors; no finite-N CAR check in this script',environment={'python':platform.python_version()},
        words=len(rows),local_graphs=sum(r['graphs'] for r in rows),records=rows,
        input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs})
    target.write_text(json.dumps(payload,indent=2)+'\n')
    print(json.dumps({k:payload[k] for k in ['status','words','local_graphs']}),flush=True)

if __name__=='__main__':main()
