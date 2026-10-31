"""Read-only integration of the five post-Bose stage audits."""
import hashlib,json
from decimal import Decimal as D,localcontext


def audit(root):
    from audit_ds3_collision import audit as collision
    from audit_ds3_full_tail import audit as full_tail
    from audit_ds3_shadow_duality import audit as shadow
    import audit_ds3_peak_limit as peak
    from audit_ds3_peak_limit_refined import machin
    from audit_ds3_transverse_phase import audit as transverse
    errors=collision(root);records=[]
    for name,fn,expected in [('ds3-full-tail',full_tail,407),('ds3-shadow-duality',shadow,365),
                              ('ds3-transverse-phase',transverse,203)]:
        result=fn(root);errors.extend(name+': '+e for e in result['errors'])
        if result['checks']!=expected:errors.append(name+': unexpected check count')
        records.append(dict(stage=name,checks=result['checks'],status=result['status']))
    old=peak.PI
    try:
        peak.PI=machin()[0];result=peak.audit(root)
    finally:peak.PI=old
    errors.extend('ds3-peak-limit: '+e for e in result['errors'])
    if result['checks']!=191:errors.append('ds3-peak-limit: unexpected count')
    records.append(dict(stage='ds3-peak-limit',checks=result['checks'],status=result['status']))
    base=root/'results/stage4/ds3-transverse-phase'
    ref=json.loads((base/'refinement.json').read_text());study=json.loads((base/'study.json').read_text())
    if ref['status']!='passed' or len(ref['rows'])!=4:errors.append('transverse refinement summary')
    with localcontext() as ctx:
        ctx.prec=130
        for row in ref['rows']:
            a,b=map(D,row['result']['value']);c,d=map(D,study['rows'][row['index']]['value'])
            diff=((a-c)**2+(b-d)**2).sqrt()
            budget=D(row['result']['tail'])+D(study['rows'][row['index']]['tail'])+D('1e-75')
            if abs(diff-D(row['difference']))>D('1e-100') or diff>budget or not row['passed']:
                errors.append('transverse refinement difference '+str(row['index']))
    for name,sha in ref['hashes'].items():
        if hashlib.sha256((root/name).read_bytes()).hexdigest()!=sha:errors.append('refinement hash '+name)
    records.append(dict(stage='ds3-collision',independent_100_digit_comparisons=12,
                        scope='saved arithmetic and recorded empirical allowances'))
    return errors,records
