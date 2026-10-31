"""Independent n=0 A.5 three-mode integral audit for arXiv:2411.08090v1.

Uses the literal published A.7 and a separately labelled counterfactual
which only divides it by Gamma(1/2+i*s3)^2. No constants are fitted.
CG is computed both from the published E.22-adjacent explicit expression
and from E.15/E.17's absolutely convergent Fourier convolution.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
from time import perf_counter
import mpmath

ROOT = Path(__file__).resolve().parents[1]


def pack(c,z):
    z=c.mpc(z)
    return {'real':c.nstr(z.real,c.dps), 'imag':c.nstr(z.imag,c.dps),
            'abs':c.nstr(abs(z),c.dps), 'arg':c.nstr(c.arg(z),c.dps)}


def psi0(c,s,t,derivatives=False):
    a=c.mpf('.5')-c.j*s; b=c.conj(a)
    z=(1-c.j*c.sinh(t))/2
    q=1/(2*c.sqrt(c.cosh(c.pi*s)))
    f=q*c.hyp2f1(a,b,1,z)
    if not derivatives:
        return f
    zp=-c.j*c.cosh(t)/2; zpp=-c.j*c.sinh(t)/2
    fz=q*a*b*c.hyp2f1(a+1,b+1,2,z)
    fzz=q*a*b*(a+1)*(b+1)/2*c.hyp2f1(a+2,b+2,3,z)
    return f, fz*zp, fzz*zp*zp+fz*zpp


def observer0(c,s,t):
    return c.sqrt(s*c.tanh(c.pi*s)/c.pi)*psi0(c,s,t)


def literal_vertex(c,s1,s2,s3):
    four=[c.mpf('.25')+c.j*(-s1+s2+s3)/2,
          c.mpf('.25')+c.j*(-s1-s2+s3)/2,
          c.mpf('.25')+c.j*(-s1+s2-s3)/2,
          c.mpf('.25')+c.j*(s1+s2+s3)/2]
    return (-c.exp(c.pi*(-5*c.j+2*(s1+s2+s3))/4)*c.cosh(c.pi*s1)
            *c.fprod(c.gamma(z) for z in four)*c.gamma(c.mpf('.5')+c.j*s3)
            *c.sqrt(c.fprod(s*c.tanh(c.pi*s) for s in [s1,s2,s3]))
            /(2*c.sqrt(2)*c.pi*c.sqrt(c.fprod(1+c.exp(2*c.pi*s) for s in [s1,s2,s3]))
              *c.gamma(c.mpf('.5')-c.j*s1)*c.gamma(c.mpf('.5')+c.j*s2)))


def cg_fourier(c,d1,d2,d3):
    # E.15: the three sine exponents; Fourier f_alpha(k) is even in k.
    alphas=[d1+d3-d2,d2+d3-d1,d1+d2-d3]
    f0=c.fprod(2**a*c.gamma(1-a)/c.gamma(1-a/2)**2 for a in alphas)
    upper=[1]+[a/2 for a in alphas]
    lower=[1-a/2 for a in alphas]
    hyper=c.hyper(upper,lower,1,strict=True,sum_method='r+s+e')
    return f0*(2*hyper-1)


def cg_published(c,d1,d2,d3):
    # n1=n2=n3=0 in the explicit C_S formula immediately after E.22.
    x=(-d1+d2+d3)/2
    y=(d1-d2+d3)/2
    z=(d1+d2-d3)/2
    def reg(upper,lower):
        return c.hyper(upper,lower,1,strict=True,sum_method='r+s+e')/c.fprod(c.gamma(v) for v in lower)
    term1=((c.cos(c.pi*d3)-c.cos(c.pi*(d2-d1)))
           *c.gamma(x+1)*c.gamma(y+1)/c.gamma(-z)
           *reg([1,z+1,y+1,x+1],[2-x,2-y,2-z]))
    term2=((c.cos(c.pi*d2)-c.cos(c.pi*(d3-d1)))
           *c.gamma(x)*c.gamma(z)/c.gamma(1-y)
           *reg([1,z,y,x],[1-x,1-y,1-z]))
    pref=(2**(d1+d2+d3-1)*c.gamma(1-2*x)*c.gamma(1-2*y)*c.gamma(1-2*z)/c.pi**2)
    return pref*(term1-term2)


def tail_bound(c,labels,tmax):
    # 2F1 connection at infinity; coefficient ratio in its 1/z series
    # is <= B_s=max(1,1/4+s^2). For t>=T, |1/z|<=2/cosh(T).
    T=c.mpf(tmax); theta=c.arg(-(1-c.j*c.sinh(T))/2)
    coeff=[]
    for s in labels:
        B=max(c.mpf(1),c.mpf('.25')+s*s)
        denominator=1-2*B/c.cosh(T)
        if denominator<=0:
            raise ValueError('Tail geometric envelope requires a larger T')
        D=abs(c.gamma(2*c.j*s)/c.gamma(c.mpf('.5')+c.j*s)**2)
        L=(c.sqrt(s*c.tanh(c.pi*s)/c.pi)*D*c.sqrt(2)
           *c.cosh(s*theta)/c.sqrt(c.cosh(c.pi*s))/denominator)
        coeff.append(L)
    return 8*c.sqrt(2)*c.pi*c.fprod(coeff)*c.exp(-T/2)


def integrate(c,labels,order,cutoffs,panel_width=2):
    nodes,weights=c.gauss_quadrature(order,'legendre')
    running=c.mpc(0); result=[]; cutoff_set=set(cutoffs)
    for lo in range(0,max(cutoffs),panel_width):
        hi=lo+panel_width
        subtotal=c.mpc(0)
        for x,w in zip(nodes,weights):
            t=(lo+hi)/2+c.mpf(panel_width)*x/2
            values=[observer0(c,s,t) for s in labels]
            subtotal += w*c.cosh(t)*c.conj(values[0])*values[1]*values[2]
        running += c.mpf(panel_width)*subtotal/2
        if hi in cutoff_set:
            value=4*c.pi*c.re(running)
            result.append({'tmax':hi,'integral':pack(c,value),
                           'positive_half_integral':pack(c,running),
                           'analytic_tail_bound_evaluated':c.nstr(tail_bound(c,labels,hi),c.dps)})
    return result


def one_run(label_text,dps,order,cutoffs):
    start=perf_counter(); c=mpmath.mp.clone(); c.dps=dps
    labels=list(map(c.mpf,label_text)); s1,s2,s3=labels
    # E.19: inputs s3,s2, output s1, then E.24 n=0 phase exp(-i*pi/4).
    dims=[c.mpf('.5')-c.j*s3,c.mpf('.5')-c.j*s2,c.mpf('.5')+c.j*s1]
    cf=cg_fourier(c,*dims); cp=cg_published(c,*dims)
    phase=c.exp(-c.j*c.pi/4)
    cg=phase*cf
    literal=literal_vertex(c,*labels)
    counterfactual=literal/c.gamma(c.mpf('.5')+c.j*s3)**2
    checks=[]
    for s in labels:
        for t in map(c.mpf,['0','.7','2']):
            f,fp,fpp=psi0(c,s,t,True)
            kg=2*c.pi*c.j*c.cosh(t)*(c.conj(f)*fp-f*c.conj(fp))
            ode=fpp+c.tanh(t)*fp+(c.mpf('.25')+s*s)*f
            checks.append({'s':str(s),'t':str(t),'kg_norm':pack(c,kg),
                           'kg_error':c.nstr(abs(kg-1),dps),'ode_error':c.nstr(abs(ode),dps),
                           'time_reflection_error':c.nstr(abs(psi0(c,s,-t)-c.conj(f)),dps)})
    vals=integrate(c,labels,order,cutoffs)
    measured=c.mpf(vals[-1]['integral']['real'])
    pred_l=literal*c.conj(cg)
    pred_cf=counterfactual*c.conj(cg)
    return {'dps':dps,'gauss_order_per_width_2_panel':order,'labels':label_text,
            'mode_checks':checks,'cg_fourier_raw':pack(c,cf),'cg_published_raw':pack(c,cp),
            'cg_relative_difference':c.nstr(abs(cf-cp)/abs(cf),dps),
            'cg_input_swap_error':c.nstr(abs(cf-cg_fourier(c,dims[1],dims[0],dims[2]))/abs(cf),dps),
            'cg_E24':pack(c,cg),'literal_A7':pack(c,literal),'counterfactual_A7':pack(c,counterfactual),
            'literal_D8_prediction':pack(c,pred_l),'counterfactual_D8_prediction':pack(c,pred_cf),
            'cutoff_study':vals,'literal_relative_discrepancy':c.nstr(abs(measured-pred_l)/abs(measured),dps),
            'counterfactual_relative_discrepancy':c.nstr(abs(measured-pred_cf)/abs(measured),dps),
            'literal_prediction_over_integral':pack(c,pred_l/measured),
            'counterfactual_prediction_over_integral':pack(c,pred_cf/measured),
            'elapsed_seconds':perf_counter()-start}


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--quick',action='store_true'); args=ap.parse_args()
    output=ROOT/'results/stage3/cg-domain-source-audit'/('vertex-check-pilot.json' if args.quick else 'vertex-check.json')
    configs=[(25,20,[16,24,32])] if args.quick else [(30,24,[16,32,48,64,80]),(40,36,[80])]
    result={'scope':'n1=n2=n3=0; three distinct real principal-series labels; no fitted constants',
            'source':'https://arxiv.org/html/2411.08090v1',
            'source_equations':['A.5','A.6','4.9','D.8','E.15','E.17','E.19','E.22-adjacent C_S','E.24','A.7'],
            'tail_status':'Analytic inequality evaluated with high precision; not outward-rounded interval arithmetic',
            'quadrature_status':'Order/precision refinement is empirical, not a certified quadrature enclosure',
            'quick':args.quick,'runs':[]}
    for dps,order,cutoffs in configs:
        run=one_run(['0.4','0.7','1.2'],dps,order,cutoffs)
        result['runs'].append(run)
        result['code_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
        print(json.dumps({'dps':dps,'order':order,'seconds':run['elapsed_seconds'],
                          'cg_difference':run['cg_relative_difference'],
                          'integral':run['cutoff_study'][-1]['integral'],
                          'tail_bound':run['cutoff_study'][-1]['analytic_tail_bound_evaluated'],
                          'literal_discrepancy':run['literal_relative_discrepancy'],
                          'counterfactual_discrepancy':run['counterfactual_relative_discrepancy']},ensure_ascii=False),flush=True)


if __name__=='__main__':
    main()