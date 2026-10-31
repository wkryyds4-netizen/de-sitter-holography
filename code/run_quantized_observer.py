"""Reproducible finite preparation tasks; no edits to recoil data.

Run from any directory. Outputs results/stage3/packets/study.json and
validation.json; every tolerance is an implementation check, not a proof.
"""
from __future__ import annotations
import hashlib
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path
import mpmath as mp
import numpy as np
import scipy
from scipy.integrate import quad
from quantized_observer import (
    DoublePacket, mass, s_of_p, p_of_s, dp_ds, rho_pl, gauss_interval,
    log_gamma8, hh_kernel, right_kernel, lr_kernel, left_right_packet_matrix,
    right_two_point, right_two_point_adaptive,
)

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results'/'stage3'/'packets'
CHECKS=[]


def serial(value):
    if isinstance(value,np.ndarray): return serial(value.tolist())
    if isinstance(value,complex): return {'real':value.real,'imag':value.imag}
    if isinstance(value,np.generic): return serial(value.item())
    if isinstance(value,dict): return {str(k):serial(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)): return [serial(v) for v in value]
    return value


def check(name,passed,**details):
    CHECKS.append(dict(name=name,passed=bool(passed),**serial(details)))


def close(name,value,reference,atol=2e-11,rtol=2e-10):
    val,ref=np.asarray(value),np.asarray(reference)
    err=float(np.max(np.abs(val-ref)))
    scale=float(np.max(np.abs(ref)))
    tol=atol+rtol*scale
    check(name,err<=tol,maximum_absolute_error=err,tolerance=tol,reference_scale=scale)
    return err


def mp_gamma8(nu,a,b):
    return mp.fprod(mp.gamma(mp.mpf('0.25')+mp.j*(sn*nu+sa*a+sb*b)/2)
                    for sn in (-1,1) for sa in (-1,1) for sb in (-1,1))


def main():
    packet=DoublePacket()
    times=np.array([-2.,-1.,-.5,0.,.5,1.,2.])
    mp.mp.dps=70
    hh_integral=mp.quad(lambda s:s*mp.tanh(mp.pi*s)/(4*mp.pi**2*mp.cosh(mp.pi*s)**2),
                        [0,1,mp.inf])
    hh_exact=1/(8*mp.pi**4)
    check('HH norm: independent 70-digit integral',abs(hh_integral-hh_exact)<mp.mpf('1e-60'),
          integral=mp.nstr(hh_integral,65),exact=mp.nstr(hh_exact,65),
          absolute_error=mp.nstr(abs(hh_integral-hh_exact),8),threshold='1e-60')
    kernel_refs=[]
    for nu in (.5,1.,1.3):
        for a,b in ((.2,.3),(1.,2.),(5.,8.),(50.,53.)):
            log_ref=mp.re(mp.log(mp_gamma8(mp.mpf(str(nu)),mp.mpf(str(a)),mp.mpf(str(b)))))
            value=float(log_gamma8(nu,a,b))
            close(f'Gamma8 log vs eight direct mp factors {nu,a,b}',value,float(log_ref),
                  atol=2e-11,rtol=2e-13)
            d=float(hh_kernel(nu,a,b))
            v=float(right_kernel(nu,a,b))
            lr=float(lr_kernel(nu,a,b))
            close(f'D symmetry {nu,a,b}',d,float(hh_kernel(nu,b,a)),atol=1e-290,rtol=2e-11)
            close(f'V from D/rho {nu,a,b}',v,d*np.cosh(np.pi*b)**2/float(rho_pl(b)),
                  atol=1e-100,rtol=2e-11)
            close(f'LR from D {nu,a,b}',lr,d*np.cosh(np.pi*a)*np.cosh(np.pi*b),
                  atol=1e-100,rtol=2e-11)
            kernel_refs.append(dict(nu=nu,a=a,b=b,log_gamma8=value,D=d,V=v,LR=lr))

    p,wp=packet.nodes(96)
    pplus=packet.probability(p,1)
    pminus=packet.probability(p,-1)
    close('real signs have identical probability density',pplus,pminus,atol=0,rtol=0)
    close('p-packet normalization plus',np.dot(wp,pplus),1.)
    close('p-packet normalization minus',np.dot(wp,pminus),1.)
    close('packet energy mean',np.dot(wp,p*pplus),6.)
    for sign in (-1,1):
        dnorm=float(np.dot(wp,packet.amplitude(p,sign,True)**2))
        close(f'clock derivative norm sign {sign}',dnorm,np.pi**2/3)
    mapped_norms=[]
    for m0 in (.5,2.,5.):
        total=0.
        for lo,hi in packet.intervals:
            ss,ww=gauss_interval(float(s_of_p(lo,m0)),float(s_of_p(hi,m0)),96)
            f=packet.physical_wavefunction(ss,m0,1)
            total+=float(np.dot(ww,rho_pl(ss)*f**2))
        close(f'physical Plancherel norm after Jacobian m0={m0}',total,1.)
        mapped_norms.append(dict(m0=m0,norm=total))
    # Delta p is a variable change, not an approximation to s or m.
    close('exact mass and flat clock Hamiltonian',mass(s_of_p(p,2.)),2.+p/(2*np.pi))

    right_rows=[]
    max_refine=0.
    max_extend=0.
    for m0 in (.5,2.,5.):
        for nu in (.5,1.):
            for sigma in (.25,.5):
                coarse=right_two_point(packet,m0,nu,sigma,times,32,24,12.)
                fine=right_two_point(packet,m0,nu,sigma,times,64,40,12.)
                extended=right_two_point(packet,m0,nu,sigma,times,64,40,16.)
                vals=fine['values']
                tag=f'm0={m0},nu={nu},sigma={sigma}'
                err=close('right quadrature refinement '+tag,coarse['values'],vals,
                          atol=2e-11,rtol=2e-10)
                tail=close('right cutoff extension '+tag,extended['values'],vals,
                           atol=2e-11,rtol=2e-10)
                max_refine=max(max_refine,err)
                max_extend=max(max_extend,tail)
                close('right conjugate symmetry '+tag,vals,np.conj(vals[::-1]))
                c0=float(vals[3].real)
                check('right positive C0 '+tag,c0>0 and abs(vals[3].imag)<1e-13,C0=c0)
                check('right positive spectral magnitude bound '+tag,
                      np.max(np.abs(vals))<=c0+2e-11,
                      maximum=float(np.max(np.abs(vals))),C0=c0)
                right_rows.append(dict(m0=m0,nu=nu,sigma_proper=sigma,times_proper=times,
                    values=vals,normalized=vals/c0,coarse_values=coarse['values'],
                    maximum_refinement_difference=err,maximum_cutoff_extension_difference=tail,
                    a_cutoff=fine['a_cutoff'],extended_a_cutoff=extended['a_cutoff'],
                    p_nodes=fine['p_nodes'],a_nodes=fine['a_nodes'],
                    state_signs=[1,-1],same_values_reason='External probability is exactly identical; (5.8) is diagonal.'))
                print('right',tag,'C0=',format(c0,'.12g'),flush=True)

    left_rows=[]
    max_lr=0.
    for m0 in (.5,2.,5.):
        for nu in (.5,1.):
            coarse=left_right_packet_matrix(packet,m0,nu,32)
            fine=left_right_packet_matrix(packet,m0,nu,64)
            tag=f'm0={m0},nu={nu}'
            err=close('LR quadrature refinement '+tag,coarse['lobe_matrix'],fine['lobe_matrix'],
                      atol=2e-12,rtol=2e-10)
            max_lr=max(max_lr,err)
            close('LR real symmetric lobe matrix '+tag,fine['lobe_matrix'],fine['lobe_matrix'].T)
            close('LR real-sign interference identity '+tag,fine['contrast'],2*fine['lobe_matrix'][0,1])
            check('LR nonzero positive sign contrast '+tag,fine['contrast']>1e-10,
                  contrast=fine['contrast'])
            left_rows.append(dict(m0=m0,nu=nu,**fine,maximum_refinement_difference=err,
                                  interpretation='two-observer propagation; no positive-operator assumption'))
            print('LR',tag,'plus/minus=',fine['expectations'],flush=True)

    print('independent adaptive right two-point beginning',flush=True)
    adaptive=right_two_point_adaptive(packet,2.,1.,.5,[0.,1.])
    reference=right_two_point(packet,2.,1.,.5,[0.,1.],64,40,16.)
    adaptive_diff=close('independent nested adaptive quadrature to infinite a',adaptive['values'],
                       reference['values'],atol=2e-9,rtol=2e-8)
    adaptive_record=dict(m0=2.,nu=1.,sigma_proper=.5,times_proper=[0.,1.],**adaptive,
                         production_values=reference['values'],maximum_difference=adaptive_diff,
                         certification='adaptive error estimates and cutoff-extension checks are not interval bounds')

    # A source-audit diagnostic: pointwise positivity does not imply a PSD LR operator.
    anu=mp.mpf('1.3'); xa=mp.mpf('.1'); xb=mp.mpf('1.3')
    def mp_lr(a,b):
        return mp.re(a*mp.sinh(mp.pi*a)*b*mp.sinh(mp.pi*b)*mp_gamma8(anu,a,b)/(64*mp.pi**7))
    aa,ab,bb=mp_lr(xa,xa),mp_lr(xa,xb),mp_lr(xb,xb)
    det=aa*bb-ab*ab
    eigmin=(aa+bb-mp.sqrt((aa-bb)**2+4*ab*ab))/2
    check('LR point-kernel indefinite diagnostic at 70 digits',det<0 and eigmin<0,
          determinant=mp.nstr(det,45),smallest_eigenvalue=mp.nstr(eigmin,45))
    lr_indefinite=dict(nu='1.3',points_s=['0.1','1.3'],matrix=[[mp.nstr(aa,45),mp.nstr(ab,45)],
                      [mp.nstr(ab,45),mp.nstr(bb,45)]],determinant=mp.nstr(det,45),
                      smallest_eigenvalue=mp.nstr(eigmin,45),
                      caveat='high-precision numerical point-kernel diagnostic; not an interval proof or the chosen packet pair')

    files=[Path(__file__),Path(__file__).with_name('quantized_observer.py')]
    hashes={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in files}
    metadata=dict(created_utc=datetime.now(timezone.utc).isoformat(),python=platform.python_version(),
                  numpy=np.__version__,scipy=scipy.__version__,mpmath=mp.__version__,code_sha256=hashes,
                  source='https://arxiv.org/html/2411.08090v1',equations=['4.31','4.35','5.8','5.6','6.4','6.17'],
                  model='dS2 principal-series free matter, fully quantized observers, strict G->0, L=1',
                  gauging='SO+(2,1); selected real packets are also T invariant',
                  mass_cutoff='Q_m0 at each right-field insertion; intermediate s>=sqrt(m0^2-1/4)',
                  units='sigma and tau are proper time in L units; p=2*pi*(m-m0)',
                  kernel_scope='right V is a spectral weight; LR is not assumed positive',
                  tolerances='recorded per check; no rigorous interval quadrature claim')
    study=dict(metadata=metadata,packet=dict(c1=3.,c2=9.,width=1.,signs=[1,-1],
                  mean_p=6.,canonical_clock_time_rms=float(np.pi/np.sqrt(3)),
                  identical_resources='same |zeta|^2 and same integral |zeta prime|^2'),
               mapped_norms=mapped_norms,hh_norm=mp.nstr(hh_integral,65),kernel_references=kernel_refs,
               right_two_point=right_rows,left_right_propagation=left_rows,
               independent_quadrature=adaptive_record,lr_indefinite_diagnostic=lr_indefinite)
    failures=[c for c in CHECKS if not c['passed']]
    validation=dict(metadata=metadata,total=len(CHECKS),passed=len(CHECKS)-len(failures),failed=len(failures),
                    maxima=dict(right_refinement=max_refine,right_cutoff_extension=max_extend,
                                lr_refinement=max_lr,independent_adaptive_difference=adaptive_diff),checks=CHECKS)
    OUT.mkdir(parents=True,exist_ok=True)
    for filename,data in [('study.json',study),('validation.json',validation)]:
        (OUT/filename).write_text(json.dumps(serial(data),ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({'total':len(CHECKS),'passed':len(CHECKS)-len(failures),'failed':len(failures),
                      'maxima':validation['maxima'],'output':str(OUT)},indent=2),flush=True)
    if failures:
        print(json.dumps(failures,indent=2))
        raise SystemExit(1)

if __name__=='__main__':
    main()