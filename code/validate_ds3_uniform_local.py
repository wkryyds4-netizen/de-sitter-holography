"""Separate unmerged representation at higher precision, including tiny offsets."""
import hashlib,json
from pathlib import Path
import mpmath as mp
from ds3_shadow_duality import physical
from ds3_uniform_local import regularized_logs
from run_ds3_collision import clean


def main():
    root=Path(__file__).resolve().parents[1];base=root/'results/stage4/ds3-uniform-local';dest=base/'validation.json'
    if dest.exists():raise RuntimeError('Preserve prior validation')
    study=json.loads((base/'study.json').read_text());rows=[];mp.mp.dps=140
    for index in [0,2,3,4,6,9,10,11,16,17,18]:
        old=study['rows'][index];M,u,v,p,q,r=map(mp.mpf,[old[k] for k in ['M','u','v','p','q','r']])
        value,oldtail,detail=physical(str(M+u),str(M+v+r),[str(M+u+v),str(q),str(M),str(p)],4)
        diff=abs(value-mp.mpc(*old['value']));co=old['Ydetails']['constants']
        N=detail['omega'][1]['details']['groups'][0]['N'];m=mp.mpf(co['m']);C0=mp.mpf(co['C0']);Q=mp.mpf(co['Q'])
        product=mp.fprod(mp.mpf(j*j)/(j*j+m*m) for j in range(1,N+1))
        ey=2*C0*product/N*(Q+7*(2+mp.log(N)))
        X,Y=[o['value'] for o in detail['omega']];ex=detail['omega'][0]['tail']
        improved_tail=sum(map(abs,detail['coefficients']))*(abs(X)*ey+abs(Y)*ex+ex*ey)
        la,lb,A,z,b,d0=regularized_logs(M,u,v,p,q,r)
        aa,bb=mp.exp(la),mp.exp(lb)
        da=aa*.5j*(-mp.digamma(1+z)-mp.digamma(1-1j*b-z)+sum(mp.digamma(a-z) for a in A[2:]))
        db=bb*.5j*(mp.digamma(1-z)+2*mp.digamma(1+1j*(d0+r))-mp.digamma(1+1j*d0-z)-sum(mp.digamma(a-z) for a in A[:2]))
        row=dict(index=index,reference=value,reference_tail=improved_tail,old_unmerged_tail=oldtail,
                 difference=diff,relative=diff/abs(value),rounding_allowance=mp.mpf('1e-65'),
                 passed=diff<=mp.mpf(old['tail'])+improved_tail+mp.mpf('1e-65'),
                 regularized_coefficients=[aa,bb],coefficient_derivatives=[da,db],
                 coefficient_bound_passed=max(abs(aa),abs(bb))<=C0,
                 derivative_bound_passed=abs(da)+abs(db)<=C0*Q,reference_N=N)
        rows.append(clean(row));print(index,mp.nstr(diff,5),row['passed'],flush=True)
    paths=['code/validate_ds3_uniform_local.py','code/ds3_uniform_local.py','code/ds3_shadow_duality.py','results/stage4/ds3-uniform-local/study.json']
    result=dict(status='passed' if all(r['passed'] and r['coefficient_bound_passed'] and r['derivative_bound_passed'] for r in rows) else 'failed',
                rows=rows,dps=140,scope='Independent unmerged sums at doubled cutoff; near-zero analytic merged tail applied to identical truncations. Rounding allowance empirical.',
                hashes={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in paths})
    dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print(result['status'])


if __name__=='__main__':main()
