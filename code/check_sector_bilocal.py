"""Fresh sector-interface evidence; refuse to overwrite a successful output."""
from datetime import datetime,timezone
from pathlib import Path
import hashlib,json,platform,time
import numpy as np
import scipy
import sector_bilocal as sb

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/susskind-sector-bilocal/interface-checks.json'

def main():
    if OUT.exists():raise FileExistsError(OUT)
    start=time.time();colored=sb.colored_checks()
    print(json.dumps({'colored':'passed'}),flush=True)
    number=sb.number_rows()
    charged=[sb.charged_rows(*args) for args in [(3,1,2,3,811),(4,1,2,2,812),(4,2,3,2,813)]]
    print(json.dumps({'charged':'passed'}),flush=True)
    bilocal=[sb.bilocal_rows(kind) for kind in ['CAR','degenerate']]
    counts=dict(source_pairings=len(colored['source_rows']),
                colored_words=sum(x['colored_words'] for x in colored['parameters']),
                replica_words=sum(x['whole_block_replica'] for x in colored['parameters']),
                number_remainders=len(number),charged_rows=sum(len(x['rows']) for x in charged),
                bilocal_times=sum(len(x['rows']) for x in bilocal),
                bilocal_bounds=3*sum(len(x['rows']) for x in bilocal))
    inputs=['code/check_sector_bilocal.py','code/sector_bilocal.py','code/sector_memory_return.py',
            'code/complex_syk_dictionary.py','code/slater_joint_site_polynomial.py',
            'code/slater_composite_sources.py','code/memory_body_source_chords.py']
    result=dict(status='passed',created_at_utc=datetime.now(timezone.utc).isoformat(),seconds=time.time()-start,
                counts=counts,colored=colored,number=number,charged=charged,bilocal=bilocal,
                parameters=dict(q=['1/4','2/3'],C_m=1,quadrature_orders=[10,18],
                                interpretation='Exact rational algebra plus finite float64 identities; not a certified quadrature or asymptotic fit.'),
                environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
                input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs})
    OUT.parent.mkdir(parents=True,exist_ok=True)
    with OUT.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(dict(status='passed',counts=counts,sha256=hashlib.sha256(OUT.read_bytes()).hexdigest())))

if __name__=='__main__':main()
