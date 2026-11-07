"""Read the external ZIP as data; independently check rational certificate assembly.

No external script execution, no extraction, no network, no modification of old evidence.
This does NOT reintegrate the 21,143 cells or prove the cubature estimates.
"""
from pathlib import Path, PurePosixPath
from fractions import Fraction as Q
from math import isqrt
from collections import Counter
import zipfile,json,hashlib,io,platform
import numpy as np

ROOT=Path(__file__).resolve().parent.parent
SRC=ROOT/'sources/external-reviews/issue-1'
OUT=ROOT/'results/issue-1-review'
def sha(b): return hashlib.sha256(b).hexdigest()
def q(x): return Q(*map(int,x)) if isinstance(x,list) else Q(x)
def end(d,k): return q(d[k])
def norm(a): return sum((Q(float(z.real))**2+Q(float(z.imag))**2 for z in a.flat),Q(0))
def sqrt_bounds(x,bits=320):
    assert x>=0
    scale=1<<bits
    k=isqrt((x.numerator*scale*scale)//x.denominator)
    lo,hi=Q(k,scale),Q(k+1,scale)
    assert lo*lo<=x<=hi*hi
    return lo,hi
def enc(x): return dict(exact=str(x),display=format(float(x),'.17g'))

if __name__=='__main__':
    dest=OUT/'independent-assembly.json'
    if dest.exists(): raise FileExistsError(dest)
    z=zipfile.ZipFile(SRC/'pentagon-review-evidence.zip'); prefix='pentagon-review-evidence/'
    assert z.testzip() is None
    names=z.namelist(); assert len(names)==len(set(names))
    assert all(not PurePosixPath(n).is_absolute() and '..' not in PurePosixPath(n).parts for n in names)
    read=lambda n:z.read(prefix+n)
    obj=lambda n:json.loads(read(n))
    manifest=obj('PACKAGE-SHA256.json')
    for n,digest in manifest.items(): assert sha(read(n))==digest,n
    assert set(names)=={prefix+'PACKAGE-SHA256.json'}|{prefix+n for n in manifest}
    assert sha((SRC/'pentagon-query.pdf').read_bytes())==sha(read('pentagon-query.pdf'))
    upstream=obj('upstream-files.json')
    for n,digest in upstream['sha256'].items(): assert sha((ROOT/n).read_bytes())==digest,n
    with np.load(io.BytesIO(read('upstream/results/stage4/flip-domain/coefficients-06.npz')),allow_pickle=False) as d:
        left,right=d['left'],d['right']
    for a in [left,right]: assert a[0,0,0]==1 and np.count_nonzero(a[0])==1 and np.isfinite(a).all()
    records=[]; total=Q(0); recomposed_total=Q(0)
    for name,branch,step in [('tensor-gaussian','right',1),('tensor-left-1','left',1),('tensor-left-2','left',2),('tensor-right-3','right',3),('tensor-right-5','right',5)]:
        fresh=obj('fresh-'+name+'.json'); arr=left if branch=='left' else right
        nv,nw=norm(arr[step-1]),norm(arr[step])
        overlap_lower=end(fresh['sum_real'],'lower')-sum((end(fresh[k],'upper') for k in ['cubature_error','space_tail','product_tail']),Q(0))
        e2=nv+nw-2*overlap_lower
        new_upper=sqrt_bounds(e2)[1]; saved=end(fresh['residual_upper'],'upper')
        # Stored residual need not enclose independently rounded endpoint arithmetic.
        # Recompute and use our own larger conservative bound in every case.
        total+=saved; recomposed_total+=new_upper
        rows=[json.loads(line) for line in read('upstream/results/stage4/residual-certificate/'+name+'-cells.jsonl').splitlines()]
        assert len(rows)==fresh['cells']; assert fresh['disjoint_saved_cells']==[]
        records.append(dict(name=name,cells=len(rows),squared_error_upper=enc(e2),
                            saved_upper=enc(saved),independent_upper=enc(new_upper),
                            saved_upper_squared_covers_independent_e2=saved*saved>=e2))
    chirp=obj('replayed-chirp.json')
    assert chirp['input_sha256']==sha((ROOT/chirp['input_file']).read_bytes())
    for row in chirp['chirp_steps']:
        step=row['step']; nv,nw=norm(right[step-1]),norm(right[step])
        e2=nv+nw-2*end(row['overlap_real'],'lower')
        new_upper=sqrt_bounds(e2)[1]; saved=end(row['residual'],'upper')
        total+=saved; recomposed_total+=new_upper
        records.append(dict(name='chirp-'+str(step),squared_error_upper=enc(e2),saved_upper=enc(saved),
                            independent_upper=enc(new_upper),saved_upper_squared_covers_independent_e2=saved*saved>=e2))
    a,b=left[-1],right[-1]; na,nb=norm(a),norm(b)
    real=Q(0); imag=Q(0)
    for x,y in zip(a.flat,b.flat):
        xr,xi,yr,yi=map(Q,[float(x.real),float(x.imag),float(y.real),float(y.imag)])
        real+=yr*xr+yi*xi; imag+=yr*xi-yi*xr
    sep2=na+nb-2*real; sep_lo,sep_hi=sqrt_bounds(sep2)
    ip_lo,ip_hi=sqrt_bounds(real*real+imag*imag)
    phase2_lo=na+nb-2*ip_hi; phase2_hi=na+nb-2*ip_lo
    phase_lo=sqrt_bounds(phase2_lo)[0];phase_hi=sqrt_bounds(phase2_hi)[1]
    margin=sep_lo-recomposed_total; phase_margin=phase_lo-recomposed_total
    assert margin>Q('0.04489650872')
    assert phase_margin>Q('0.04489650465')
    summary=obj('fresh-certificate-summary.json')
    checks={}
    checks['saved_separation_encloses_exact_squared_distance']=end(summary['separation'],'lower')**2<=sep2<=end(summary['separation'],'upper')**2
    checks['saved_total_encloses_exact_sum']=end(summary['total_step_error'],'lower')<=total<=end(summary['total_step_error'],'upper')
    checks['saved_phase_separation_intersects_independent_enclosure']=end(summary['phase_optimized_finite_separation'],'lower')<=phase_hi and phase_lo<=end(summary['phase_optimized_finite_separation'],'upper')
    assert all(checks.values())
    audit=obj('replayed-audit.json'); failures=Counter(s.split(' bound ')[0] if ' bound ' in s else 'OTHER' for s in audit['failed'])
    assert len(audit['failed'])==9890 and failures.get('OTHER',0)==0 and audit['passed'] is False
    portability=obj('bound-portability.json')
    for row in portability['records']:
        assert row['changed_upper_endpoints']==failures[row['name']]
        assert row['changed_cell_kind']==0
    old=json.loads((ROOT/'results/stage4/residual-certificate/tensor-chain-audit.json').read_text())
    result=dict(status='passed: integrity and independent rational assembly; cell integration not rerun',
        external_issue='https://github.com/wkryyds4-netizen/de-sitter-holography/issues/1',
        package_files_verified=len(manifest),upstream_files_matching_local=len(upstream['sha256']),
        fresh_cells_reported=sum(r.get('cells',0) for r in records),
        recomputed_step_bounds=records,independent_total_error_upper=enc(recomposed_total),
        independent_norm_lower_bound=enc(margin),independent_phase_optimized_lower_bound=enc(phase_margin),
        reported_audit_failure_categories=dict(failures),summary_checks=checks,
        original_local_lower_bound=enc(end(old['word_norm_lower_bound'],'lower')),
        no_external_scripts_executed=True,full_cell_reintegration_performed=False,
        mathematical_dependency='Conditional on correctness of supplied overlap enclosures and the existing quadrature/tail estimates. Rational recomposition is independently checked here.',
        environment=dict(python=platform.python_version(),numpy=np.__version__,sqrt_enclosure_bits=320),
        input_hashes={p.relative_to(ROOT).as_posix():sha(p.read_bytes()) for p in [SRC/'issue.json',SRC/'pentagon-query.pdf',SRC/'pentagon-review-evidence.zip',Path(__file__)]})
    dest.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:result[k] for k in ['status','package_files_verified','upstream_files_matching_local','independent_norm_lower_bound','independent_phase_optimized_lower_bound']}))
