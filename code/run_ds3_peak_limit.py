"""Explicit finite-M peak bound, including the large-parameter companion series."""
import hashlib,json,platform,time
from pathlib import Path
from fractions import Fraction
import mpmath as mp
import sympy as sp
from ds3_peak_normalization import symbolic_certificates
from ds3_shadow_duality import physical
from ds3_collision_asymptotic import leading_bound
from run_ds3_collision import clean


def main():
    root=Path(__file__).resolve().parents[1];base=root/'results/stage4/ds3-peak-limit';base.mkdir(exist_ok=True)
    raw=base/'runs.jsonl';dest=base/'study.json'
    if raw.exists() or dest.exists():raise RuntimeError('Preserve prior run')
    mp.mp.dps=100;rows=[]
    for labels in [('.3','-.2','.6','.6'),('.2','.1','.7','.8')]:
        ur,vr,pr,qr=map(sp.Rational,labels);u,v,p,q=map(mp.mpf,labels)
        for M in [8,16,32,64,128,256]:
            t0=time.monotonic();m=mp.mpf(M);b=m+u;d=m+v
            value,tail,detail=physical(str(M+ur),str(M+vr),[str(M+ur+vr),str(qr),str(M),str(pr)],2)
            C=detail['coefficients'][0];D=detail['coefficients'][1]
            first=detail['omega'][0]['details'];ca,cb=first['coefficients'];D0=C*ca;E0=C*cb
            q0,RY,leading=leading_bound(M,*labels)
            Amax=max(abs(mp.im(x)) for x in leading['top'])
            coefs=[Fraction(x) for x in first['groups'][1]['certificate']['coefficients'][1:]]
            maximum=max(abs(x) for x in coefs)
            Nstar=2+(maximum.numerator+maximum.denominator-1)//maximum.denominator
            B=mp.sqrt(mp.sinh(mp.pi*b)/(mp.pi*b))*mp.exp(mp.pi*(2*b*Amax+Amax*Amax)/b)
            big_bound=(3*Nstar+3)*B
            RF=(mp.zeta(2)-1)/(1+min(b,d)**2)
            phase=mp.arg(4*mp.pi*D0);phase_difference=abs(4*mp.pi*D0-1)
            error_bound=RY/(2*mp.pi)+((RF+phase_difference)/(2*mp.pi)+2*abs(E0)*big_bound)*(abs(q0)+RY)
            leading_value=mp.re(q0)/(2*mp.pi)
            constant=(2*mp.digamma(1)-mp.re(sum(mp.digamma(x) for x in leading['top'])))/(2*mp.pi)
            expansion=mp.log(m)/mp.pi+constant+(u+v)/(2*mp.pi*m)
            theta2=(-p*p*v+q*q*u+u*u*v-u*v*v+u-v)/4
            row=clean(dict(M=M,u=u,v=v,p=p,q=q,value=value,numerical_tail=tail,leading=leading_value,leading_error=abs(value-leading_value),
                           leading_bound=error_bound,passed=abs(value-leading_value)<=error_bound+tail,
                           expansion=expansion,expansion_error=abs(value-expansion),scaled_expansion_error=abs(value-expansion)*m*m/mp.log(m),
                           constant=constant,phase=phase,phase_coefficient=theta2,phase_scaled=phase*m*m,
                           unit_modulus_residual=abs(abs(4*mp.pi*D0)-1),coefficient_conjugacy_residual=abs(D-mp.conj(C)),
                           q0=q0,RY=RY,RF=RF,Amax=Amax,Nstar=Nstar,large_series_bound=big_bound,
                           suppressed_coefficient=E0,suppressed_term_bound=abs(E0)*big_bound,leading_coefficient=D0,
                           details=detail,seconds=time.monotonic()-t0))
            with raw.open('a',encoding='utf-8') as f:f.write(json.dumps(row)+'\n');f.flush()
            rows.append(row);print(labels,M,'error',mp.nstr(abs(value-leading_value),5),'bound',mp.nstr(error_bound,5),'phase M2',mp.nstr(phase*m*m,8),flush=True)
    symbolic=symbolic_certificates()
    paths=['code/ds3_peak_normalization.py','code/run_ds3_peak_limit.py','code/ds3_shadow_duality.py','code/ds3_full_tail_bound.py','code/ds3_collision_asymptotic.py','code/ds3_local_vertex.py','code/global_racah.py']
    data=dict(status='passed' if all(r['passed'] for r in rows) and all(symbolic[k] for k in ['coefficient_conjugate_pairs','small_phase_cancels','phase_expansion_correct']) else 'failed',rows=rows,symbolic=symbolic,
              versions=dict(python=platform.python_version(),mpmath=mp.__version__,sympy=sp.__version__),dps=100,
              scope='Complete scalar kernel on the exact collision surface, fixed transfers only; excludes rounding and all mass-integral limits.',
              hashes={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in paths})
    dest.write_text(json.dumps(data,indent=2),encoding='utf-8');print(data['status'],len(rows))


if __name__=='__main__':main()
