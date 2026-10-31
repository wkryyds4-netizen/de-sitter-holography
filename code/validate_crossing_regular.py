"""Validate cyclic charts and pole-subtracted double diagonals; no large scan."""
from __future__ import annotations
from datetime import datetime,timezone
import argparse,hashlib,json,platform,sys,traceback
from pathlib import Path
import mpmath
from observer_crossing import CrossingEvaluator
from observer_crossing_regular import RegularCrossingEvaluator

ROOT=Path(__file__).resolve().parents[1]
CASES=[
 ('old_1','0.6',['0.35','0.7','1.1','1.45']),
 ('old_2','1.0',['0.4','1.25','0.8','1.7']),
 ('single_ac','0.6',['0.6','0.8','0.6','1.1']),
 ('single_bd','0.6',['0.6','0.8','1.1','0.8']),
 ('double','0.6',['0.6','0.8','0.6','0.8']),
 ('all_equal','0.6',['0.6','0.6','0.6','0.6']),
 ('double_2','1.0',['0.4','1.2','0.4','1.2']),
]


def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--output',type=Path,default=ROOT/'results/stage3/fourpoint/crossing-regular-validation.json')
 args=parser.parse_args()
 lo=RegularCrossingEvaluator(25,15)
 hi=RegularCrossingEvaluator(40,15)
 c=hi.ctx
 tol=c.mpf('1e-20')
 checks=[];rows=[];failures=[];approaches=[]
 old=json.loads((ROOT/'results/stage3/fourpoint/crossing-validation.json').read_text(encoding='utf-8-sig'))

 def rel(a,b):
  a,b=c.mpc(a),c.mpc(b)
  return abs(a-b)/max(abs(a),abs(b),c.mpf('1e-1000'))
 def complex_data(z):
  if z is None:return None
  return {'real':c.nstr(c.re(z),55),'imag':c.nstr(c.im(z),55)}
 def from_data(z):return c.mpc(z['real'],z['imag'])
 def check(name,error,limit=tol):
  error=c.mpf(error)
  checks.append(dict(name=name,passed=bool(c.isfinite(error) and error<=limit),
                     residual=c.nstr(error,30),limit=c.nstr(limit,15)))
 def flag(name,value,**details):checks.append(dict(name=name,passed=bool(value),**details))
 def pack(r):
  return dict(labels=[c.nstr(v,80) for v in r.labels],nu=c.nstr(r.matter_nu,50),
              chart_labels=[c.nstr(v,80) for v in r.chart_labels],chart_rotation=r.chart_rotation,
              method=r.method,symmetric_s=complex_data(r.symmetric_s),hh_cross=complex_data(r.hh_cross),
              raw_s_chart=complex_data(r.raw_s_chart),gamma_product_chart=complex_data(r.gamma_product_chart),
              cancellation_ratio=c.nstr(r.cancellation_ratio,30),
              residue_relative_error=c.nstr(r.residue_relative_error,30),
              diagonal_a0=complex_data(r.diagonal_a0),diagonal_b0=complex_data(r.diagonal_b0),
              diagonal_difference_derivative=complex_data(r.diagonal_difference_derivative),
              requested_dps=r.requested_dps,working_dps=r.working_dps,
              max_observed_formula_dps=r.max_observed_formula_dps,gap_guard_digits=r.gap_guard_digits,
              refinement_count=r.refinement_count,hyper_calls=r.hyper_calls,
              max_hyper_excess_error=c.nstr(r.max_hyper_excess_error,20),
              reality_relative_error=c.nstr(r.reality_relative_error,30),elapsed_seconds=r.elapsed_seconds)
 base_results={}
 for name,nu,labels in CASES:
  print('regular case '+name,flush=True)
  row=dict(name=name,nu=nu,labels=labels,orbit=[])
  try:
   base=lo.evaluate(nu,*labels);fine=hi.evaluate(nu,*labels)
   base_results[name]=fine
   row['base_low']=pack(base);row['base_high']=pack(fine)
   check(name+'/precision_S',rel(base.symmetric_s,fine.symmetric_s))
   check(name+'/precision_HH',rel(base.hh_cross,fine.hh_cross))
   check(name+'/reality',fine.reality_relative_error)
   check(name+'/residue',fine.residue_relative_error)
   if name.startswith('old_'):
    frozen=old['points'][int(name[-1])-1]['base_high']
    check(name+'/frozen_76_S',rel(fine.symmetric_s,from_data(frozen['symmetric_s'])))
    check(name+'/frozen_76_HH_constant',rel(fine.hh_cross,from_data(frozen['hh_cross'])))
    # Independently check the analytic extraction against the unstripped B.21 sum.
    with c.workdps(55):
     r=tuple(c.mpf(v) for v in labels);v=c.mpf(nu);delta=r[0]-r[2]
     aa,bb=hi.stripped_branches(r,v)
     factored=c.gamma(c.j*delta/2)*aa+c.gamma(-c.j*delta/2)*bb
     plain=c.fsum(hi.s_terms(r,v))
     check(name+'/pole_factorization_identity',rel(factored,plain),c.mpf('1e-40'))
   seen={tuple(labels):base}
   for reflected in (False,True):
    orient=list(reversed(labels)) if reflected else labels
    for shift in range(4):
     perm=orient[shift:]+orient[:shift];key=tuple(perm)
     if key not in seen:seen[key]=lo.evaluate(nu,*perm)
     result=seen[key]
     tag=f'{name}/D4/reflection_{int(reflected)}/rotation_{shift}'
     check(tag+'/symmetric_S',rel(result.symmetric_s,base.symmetric_s))
     a,b,z,d=[c.mpf(v) for v in labels];p,q,r,s=[c.mpf(v) for v in perm]
     ratio=c.cosh(c.pi*a)*c.cosh(c.pi*z)/(c.cosh(c.pi*p)*c.cosh(c.pi*r))
     check(tag+'/original_HH_prefactor',rel(base.hh_cross,ratio*result.hh_cross))
     check(tag+'/raw_reality',result.reality_relative_error)
     row['orbit'].append(dict(reflected=reflected,shift=shift,result=pack(result)))
   if name in ('double','double_2','all_equal'):
    flag(name+'/analytic_diagonal_method',base.method=='diagonal_derivative')
    flag(name+'/derivative_precision_separately_recorded',base.max_observed_formula_dps>base.working_dps)
   if name=='single_ac':flag(name+'/uses_regular_cyclic_chart',base.chart_rotation==1)
   if name=='single_bd':flag(name+'/does_not_perturb_single_bd',base.chart_rotation==0)
  except Exception as exc:
   record=dict(case=name,error_type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc())
   failures.append(record);row['exception']=record
  rows.append(row)

 # Separate full-kernel approach tests. Production diagonal values NEVER use
 # these offsets; they use the derivative of the analytic stripped terms.
 directions=[('ac',[1,0,-1,0]),('bd',[0,1,0,-1]),('mixed',[1,2,-1,-2])]
 for name,nu,labels in [CASES[4],CASES[5]]:
  if name not in base_results:continue
  target=base_results[name].symmetric_s
  selected=directions if name=='double' else [directions[0],directions[2]]
  for direction,vector in selected:
   print('approach '+name+'/'+direction,flush=True)
   entry=dict(case=name,direction=direction,vector=vector,levels=[])
   try:
    averages=[]
    for htext in ['0.0002','0.0001','0.00005']:
     h=c.mpf(htext);center=[c.mpf(v) for v in labels]
     plus=[c.nstr(x+h*v,65) for x,v in zip(center,vector)]
     minus=[c.nstr(x-h*v,65) for x,v in zip(center,vector)]
     rp=lo.evaluate(nu,*plus);rm=lo.evaluate(nu,*minus)
     avg=(rp.symmetric_s+rm.symmetric_s)/2
     averages.append(avg)
     entry['levels'].append(dict(h=htext,plus=pack(rp),minus=pack(rm),
                                average=complex_data(avg),relative_error=c.nstr(rel(avg,target),25)))
    rich1=(4*averages[1]-averages[0])/3
    rich2=(4*averages[2]-averages[1])/3
    entry['richardson_coarse']=complex_data(rich1)
    entry['richardson_fine']=complex_data(rich2)
    entry['richardson_change_relative']=c.nstr(rel(rich2,rich1),25)
    check(name+'/'+direction+'/richardson_matches_derivative',rel(rich2,target),c.mpf('1e-12'))
    check(name+'/'+direction+'/richardson_refinement',rel(rich2,rich1),c.mpf('2e-11'))
    flag(name+'/'+direction+'/symmetric_approach_improves',rel(averages[2],target)<rel(averages[0],target))
   except Exception as exc:
    record=dict(case=name,direction=direction,error_type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc())
    failures.append(record);entry['exception']=record
   approaches.append(entry)

 print('near-diagonal precision and boundary rejection',flush=True)
 near=None
 try:
  tiny=lo.evaluate('0.6','0.600000000000000000000001','0.8','0.6','0.8')
  near=pack(tiny)
  flag('tiny_nonzero_gap_not_rounded_to_diagonal',tiny.method!='diagonal_derivative' and tiny.gap_guard_digits>=24)
  check('tiny_gap_agrees_with_diagonal',rel(tiny.symmetric_s,base_results['double'].symmetric_s),c.mpf('1e-20'))
  tiny_high=hi.evaluate('0.6','0.600000000000000000000001','0.8','0.6','0.8')
  near['high_precision_repeat']=pack(tiny_high)
  check('tiny_gap_precision_refinement',rel(tiny.symmetric_s,tiny_high.symmetric_s))
 except Exception as exc:
  failures.append(dict(case='tiny_gap',error_type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc()))
 try:
  lo.evaluate('0.6','0.6'+'0'*200+'1','0.8','0.6','0.8')
 except ArithmeticError as exc:
  flag('unresolved_exact_input_gap_rejected','input gap exceeds precision budget' in str(exc),message=str(exc))
 except Exception as exc:
  flag('unresolved_exact_input_gap_rejected',False,unexpected=repr(exc))
 else:
  flag('unresolved_exact_input_gap_rejected',False,unexpected='accepted as exact diagonal')
 for name,bad in [('zero_label',('0.6','0','0.8','0.6','0.8')),('negative_nu',('-0.6','0.6','0.8','0.6','0.8'))]:
  try:lo.evaluate(*bad)
  except ValueError:flag('domain/'+name,True)
  except Exception as exc:flag('domain/'+name,False,unexpected=repr(exc))
  else:flag('domain/'+name,False,unexpected='accepted')

 failed=sum(not row['passed'] for row in checks)
 code_names=['observer_crossing.py','observer_crossing_regular.py','validate_crossing_regular.py']
 inputs=['sources/kolchmeyer-liu-2411.08090v1.html','results/stage3/fourpoint/crossing-validation.json']
 sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
 report=dict(schema='observer-crossing-regular-validation-v1',generated_utc=datetime.now(timezone.utc).isoformat(),
             status='passed' if failed==0 and not failures else 'failed',checks=len(checks),failed=failed,
             evaluation_exceptions=len(failures),check_details=checks,exceptions=failures,
             source='https://arxiv.org/html/2411.08090v1',source_equations=['6.26','6.27','A.14','B.21','B.22','B.23','B.24'],
             code_sha256={n:sha(ROOT/'code'/n) for n in code_names},
             input_sha256={n:sha(ROOT/n) for n in inputs},
             numerics=dict(low_dps=25,high_dps=40,guard_digits=15,max_working_dps=180,
                           derivative='mpmath.diff central step on analytic pole-stripped difference; addprec=15 bits; callback precision recorded separately',
                           approach_h=['0.0002','0.0001','0.00005'],interval_certification=False),
             environment=dict(python=sys.version,mpmath=mpmath.__version__,platform=platform.platform()),
             retained_initial_diagnostic='results/stage3/fourpoint/crossing-regular-initial-validation.json',
             cases=rows,approaches=approaches,tiny_gap=near,
             limitations=['Strictly positive real masses and nu; no endpoint s=0 extension.',
                          'max_working_dps caps outer context; differentiation and special functions use additional internal precision.',
                          'Analytic finite-limit formula with numerical derivatives, not an interval-certified evaluation.',
                          'Precision refinement, D4 and five approach paths are finite diagnostics, not uniform global error bounds.',
                          'No high-dimensional wavepacket scan or proof of convergence of such a scan.',
                          'The original 76-check implementation and result remain unchanged.'])
 args.output.parent.mkdir(parents=True,exist_ok=True)
 args.output.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
 print(json.dumps({k:report[k] for k in ['status','checks','failed','evaluation_exceptions']}),flush=True)
 print(args.output,flush=True)
 return 0 if report['status']=='passed' else 1

if __name__=='__main__':raise SystemExit(main())
