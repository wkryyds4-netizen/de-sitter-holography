"""Independent Decimal audit of saved Bose algebra, tails, failures and hashes."""
import hashlib
import json
from decimal import Decimal as D,localcontext
from audit_ds3_crossing_dictionary import PI


def audit(root):
    base=root/'results/stage4/ds3-bose'
    read=lambda n:json.loads((base/n).read_text(encoding='utf-8'))
    errors=[]
    def req(ok,msg):
        if not ok:errors.append('Bose: '+msg)
    def f(x):return D.from_float(x) if isinstance(x,float) else D(x)
    def z(v):return tuple(map(f,v))
    def norm(v):return (v[0]*v[0]+v[1]*v[1]).sqrt()
    def close(a,b,msg,tol='1e-18'):req(abs(a-b)<=D(tol),msg)
    summary=read('summary.json');validation=read('validation.json');repair=read('precision-repair.json');nodes=read('weighted-node-checks.json')
    records=[summary,validation,repair,nodes]
    req(summary['status']=='passed' and len(summary['results'])==3,'three replays')
    req(summary['total_points']==1216480,'total raw points')
    failed=[r for r in validation['checks'] if not r['passed']]
    req(len(validation['checks'])==22 and len(failed)==1 and validation['status']=='failed','preserved initial failure')
    req(failed[0]['name']=='physical_flip_hermiticity_regular','specific failed check')
    req(repair['status']=='passed' and len(repair['rows'])==3,'selected precision repair')
    with localcontext() as ctx:
        ctx.prec=80
        delta=lambda x:1/((1+x*x).sqrt()+x)
        sh=lambda x:(x.exp()-(-x).exp())/2
        s,c0,c1,p=map(D,['.2','.4','.8','.6'])
        C0=((PI*delta(s)).exp()-1)/(2*delta(s))
        gmin=c0*(PI*delta(c1)).exp()/(1+(-2*PI*c0).exp())
        ap=p/(4*PI*sh(PI*p));nmax=((c1+C0)*ap*ap+ap/PI**3)/gmin
        for entry in summary['results']:
            raw=read(entry['label']+'.json');records.append(raw)
            req(raw['status']=='completed','completed run')
            path=root/raw['raw_grid']
            req(hashlib.sha256(path.read_bytes()).hexdigest()==raw['raw_grid_sha256'],'raw hash')
            req(raw['points']==entry['points'],'run point count')
            req(entry['max_replay_difference']<3e-17,'independent replay residual')
            L=f(raw['parameters']['L']);decay=(-PI*(L-c1)).exp()
            ta=(PI*delta(L)).exp()*decay/(16*PI**3*gmin)*(L*L+L*(3/PI-c0)+3/PI**2-c0/PI)
            tb=decay/(8*PI**4*gmin)*(L+C0+2/PI);eps=ta+tb
            close(nmax,f(raw['mass_tail']['Nmax']),'Nmax','1e-17')
            close(eps,f(raw['mass_tail']['epsilon']),'epsilon','1e-23')
            close(4*eps+8*(eps*nmax).sqrt(),f(raw['mass_tail']['D_bound']),'D tail','1e-19')
            close(eps+2*(eps*nmax).sqrt(),f(raw['mass_tail']['bose_inner_product_bound']),'overlap tail','1e-19')
            for r in raw['rows']:
                N=f(r['N_label']);E1=z(r['E1']);E2=z(r['E2']);Y=z(r['Y']);Z=z(r['Z'])
                n1=N+E1[0];n2=N+E2[0];ov=(Y[0]+Z[0],Y[1]+Z[1])
                stable=f(r['D_direct'])+z(r['D_exchange'])[0]
                close(n1,z(r['N1'])[0],'Bose N1');close(n2,z(r['N2'])[0],'Bose N2')
                close(ov[0],z(r['overlap'])[0],'Bose overlap real');close(ov[1],z(r['overlap'])[1],'Bose overlap imag')
                close(stable,z(r['D'])[0],'stable D')
                close(n1+n2-2*ov[0],f(r['expanded_D']),'expanded D')
                close(stable,f(r['expanded_D']),'stable-expanded','1e-17')
                req(stable>=0 and n1>0 and n2>0 and norm(ov)<=(n1*n2).sqrt()+D('1e-17'),'Cauchy and positivity')
                close(norm(ov)/(n1*n2).sqrt(),f(r['cauchy_ratio']),'ratio','1e-14')
            req(raw['rows'][0]['D']==[0.,0.],'exact numerical equal-time zero')
            close(f(raw['moment_direct'])+z(raw['moment_exchange'])[0],z(raw['moment_total'])[0],'moment assembly')
            req(z(raw['moment_total'])[0]>0,'positive numerical moment')
        for r in validation['checks']:
            req(r['passed']==(r['residual']<=r['tolerance']),'check flags')
        for r in validation['spherical_derivatives']:
            req(D(r['integral'])>0,'positive spherical diagonal')
            close(D(r['integral']),D(r['paired_integral']),'paired diagonal','1e-40')
        for r in validation['kernel_references']+nodes['rows']:
            val=z(r['value'] if 'value' in r else r['u']);ref=tuple(map(D,r['reference']))
            err=norm((val[0]-ref[0],val[1]-ref[1]))
            close(err,D(r['absolute_error']),'reference absolute','1e-40')
            close(err/norm(ref),D(r['relative_error']),'reference relative','1e-35')
        req(D(validation['kernel_references'][0]['relative_error'])>1,'preserve large resonance error')
        req(D(validation['kernel_references'][-1]['relative_error'])>D('.001'),'preserve far-tail error')
        req(len(nodes['rows'])==6 and all(D(r['relative_error'])<D('1e-9') for r in nodes['rows']),'selected weighted references')
        for r in repair['rows']:
            req(r['passed'] and r['hermiticity_residual']<1e-10,'repaired selected Hermiticity')
            for key in ['precision','reverse_precision']:
                if r[key]['method']=='mpmath55_75':req(D(r[key]['refinement_relative'])<D('1e-25'),'precision refinement')
    for record in records:
        for path,sha in record['hashes'].items():
            req(hashlib.sha256((root/path).read_bytes()).hexdigest()==sha,'input hash '+path)
    return errors


if __name__=='__main__':
    from pathlib import Path
    errors=audit(Path(__file__).resolve().parents[1])
    print(json.dumps(dict(status='failed' if errors else 'passed',errors=errors),indent=2))
    raise SystemExit(bool(errors))
