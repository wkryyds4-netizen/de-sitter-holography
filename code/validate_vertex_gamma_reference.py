"""Gamma-only reference for the frozen independent n=0 vertex integral.

Bros et al. arXiv:0901.4223v1 (2),(64)-(66),(124), mapped to
Kolchmeyer-Liu arXiv:2411.08090v1 E.19/E.24. No quadrature is rerun.
"""
from pathlib import Path
import hashlib
import itertools
import json
import mpmath
from validate_vertex_source import literal_vertex, pack

ROOT=Path(__file__).resolve().parents[1]
source_path=ROOT/'results/stage3/cg-domain-source-audit/vertex-check.json'
raw=json.loads(source_path.read_text(encoding='utf8'))
c=mpmath.mp.clone(); c.dps=60
s1,s2,s3=map(c.mpf,raw['runs'][-1]['labels'])
labels=[s1,s2,s3]; quarter=c.mpf('.25'); half=c.mpf('.5')
g4=c.fprod(c.gamma(quarter+c.j*v/2) for v in [-s1+s2+s3,-s1-s2+s3,-s1+s2-s3,s1+s2+s3])
D=c.gamma(half-c.j*s1)*c.gamma(half+c.j*s2)*c.gamma(half+c.j*s3)
cs=g4/(c.pi**c.mpf('1.5')*D)
cphys=c.exp(-c.j*c.pi/4)*cs
cprod=c.fprod(c.cosh(c.pi*s) for s in labels)
sprod=c.fprod(s*c.tanh(c.pi*s) for s in labels)
g8=c.fprod(c.gamma(quarter+c.j*sum(e*s for e,s in zip(signs,labels))/2) for signs in itertools.product([-1,1],repeat=3))
h=g8*cprod/(4*c.pi**5)
integral=c.cosh(c.pi*s1)*c.sqrt(sprod)*h/(2*c.sqrt(c.pi*cprod))
lit=literal_vertex(c,*labels)
cf=lit/c.gamma(half+c.j*s3)**2
pred=cf*c.conj(cphys)
run=raw['runs'][-1]
old_cg=c.mpc(run['cg_E24']['real'],run['cg_E24']['imag'])
measured=c.mpf(run['cutoff_study'][-1]['integral']['real'])
result={'source':['https://arxiv.org/html/0901.4223v1','https://arxiv.org/html/2411.08090v1','https://dlmf.nist.gov/14.24.E1'],
        'working_dps':c.dps,'labels':run['labels'],'cg_Gamma_E24':pack(c,cphys),
        'legendre_h_Gamma':pack(c,h),'n0_integral_Gamma_contour':pack(c,integral),
        'counterfactual_D8_Gamma':pack(c,pred),
        'Gamma_identity_absolute_difference':c.nstr(abs(pred-integral),c.dps),
        'original_4F3_cg_relative_difference_from_Gamma':c.nstr(abs(old_cg-cphys)/abs(cphys),c.dps),
        'frozen_integral_absolute_difference_from_Gamma':c.nstr(abs(measured-integral),c.dps),
        'frozen_integral_tail_bound':run['cutoff_study'][-1]['analytic_tail_bound_evaluated'],
        'hashes':{str(p.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),ROOT/'code/validate_vertex_source.py',source_path]}}
out=ROOT/'results/stage3/cg-domain-source-audit/vertex-analytic-reference.json'
out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
print(json.dumps(result,ensure_ascii=False,indent=2))