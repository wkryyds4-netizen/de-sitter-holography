"""Independent physical-weight reconstruction from E, W and rho; retained U grids."""
import hashlib
import json
from pathlib import Path
import numpy as np
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'results/stage4/ds3-crossed-packet'
LABELS=['pilot-n4-L8','pilot-precision-n4-L8','refined-n6-L10',
        'precision-n6-L12','refined-n8-L12']


def denom(a,b,p):
    return (np.cosh(np.pi*(a+b))+np.cosh(np.pi*p))*(np.cosh(np.pi*(a-b))+np.cosh(np.pi*p))/4


def E(a,b,p):
    return np.sqrt(np.sinh(np.pi*p)/p)/(8*np.sqrt(np.pi))*np.sinh(np.pi*a)*np.sinh(np.pi*b)*np.sqrt(np.cosh(np.pi*b)/np.cosh(np.pi*a))/denom(a,b,p)


def K(a,b,p):
    return np.pi/(8*denom(a,b,p))*np.sinh(np.pi*a)/a*np.sinh(np.pi*b)/b*np.sinh(np.pi*p)/p


def replay(label):
    data=json.loads((BASE/(label+'.json')).read_text())
    path=ROOT/data['raw_grid']
    assert hashlib.sha256(path.read_bytes()).hexdigest()==data['raw_grid_sha256']
    assert data['status']=='completed'
    with np.load(path) as f:
        x,w,c,cw,U=f['x'],f['w'],f['c'],f['cw'],f['U']
        a,b,d=np.meshgrid(x,x,x,indexing='ij');a,b,d=[v.ravel() for v in [a,b,d]]
        mass_weight=(w[:,None,None]*w[None,:,None]*w[None,None,:]).ravel()
        ns=nt=0.;values=np.zeros(4,complex);absolute=0.
        for i,(ci,wi) in enumerate(zip(c,cw)):
            ws=2*b*b/np.pi*K(a,b,.9)*K(b,ci,.6)
            wt=2*d*d/np.pi*K(ci,d,.9)*K(d,a,.6)
            # This reconstructs the F/sqrt(W) amplitudes before simplification.
            physical_weight=E(a,b,.9)*E(b,ci,.6)*E(a,d,.6)*E(d,ci,.9)/np.sqrt(ws*wt)/(2*ci*ci/np.pi)
            physical_weight*=np.exp(np.pi*(np.sqrt(1+a*a)-np.sqrt(1+ci*ci)))
            amplitude=wi*mass_weight*physical_weight*U[i]
            energy=np.sqrt(1+a*a)+np.sqrt(1+ci*ci)-np.sqrt(1+b*b)-np.sqrt(1+d*d)
            for j,t in enumerate(data['times']):values[j]+=np.sum(amplitude*np.exp(1j*t*energy))
            absolute+=float(np.sum(abs(amplitude)))
            aa=x[:,None];bb=x[None,:];ww=w[:,None]*w[None,:]
            for q,p,which in [(.9,.6,'s'),(.6,.9,'t')]:
                W=2*bb*bb/np.pi*K(aa,bb,q)*K(bb,ci,p)
                v2=E(aa,bb,q)**2*E(bb,ci,p)**2/W/(2*ci*ci/np.pi)
                n=wi*np.sum(ww*v2*np.exp(np.pi*(np.sqrt(1+aa*aa)-np.sqrt(1+ci*ci))))
                if which=='s':ns+=n
                else:nt+=n
        saved=np.array([complex(*v) for v in data['values']])
        difference=abs(values-saved)
        assert max(difference)<2e-17
        assert abs(absolute-data['absolute_sum'])<2e-17
        return dict(label=label,values=[[float(z.real),float(z.imag)] for z in values],
            replay_absolute_difference=difference.tolist(),absolute_sum=absolute,
            truncated_norms=[float(ns),float(nt)],points=int(U.size),
            initial_norm=float(sum(cw)))


def run():
    rows=[replay(label) for label in LABELS]
    data=[json.loads((BASE/(label+'.json')).read_text()) for label in LABELS]
    vals=[np.array([complex(*v) for v in d['values']]) for d in data]
    final=data[-1];ns=final['positive_norms']['ns'];nt=final['positive_norms']['nt']
    D=ns+nt-2*vals[-1].real
    bound=np.sqrt(ns*nt)
    assert np.all(abs(vals[-1])<=bound) and np.all(D>0)
    witness=json.loads((BASE/'bounds-validation.json').read_text())['witnesses'][-1]
    assert D[0]>witness['order_difference_lower_bound']
    P,Q,u,v=sp.symbols('P Q u v')
    DP=u*u+v*v-1+2*P*u*v+P*P;DQ=DP.subs(P,Q)
    numerator=2*(Q-P)*(v*u*u+(P+Q)*u+v*(1+P*Q-v*v))
    derivative_verified=sp.simplify(sp.diff(DP/DQ,u)-numerator/DQ**2)==0
    assert derivative_verified
    summary=dict(status='passed',rows=rows,final_label=LABELS[-1],
        times=final['times'],X=final['values'],Ns=ns,Nt=nt,D=D.tolist(),
        cauchy_ratios=(abs(vals[-1])/bound).tolist(),
        grid_n4_to_n6=abs(vals[2]-vals[0]).tolist(),
        grid_n6_to_n8=abs(vals[4]-vals[3]).tolist(),
        L10_to_L12_with_series_change=abs(vals[3]-vals[2]).tolist(),
        same_grid_series_change=abs(vals[1]-vals[0]).tolist(),
        witness_lower_bound=witness['order_difference_lower_bound'],
        monotonicity_symbolic_identity=bool(derivative_verified),
        total_kernel_points=sum(r['points'] for r in rows),
        empirical_only='Grid/series differences are not certified bounds on the total error.')
    paths=['code/replay_ds3_crossed_packet.py']
    paths += [f'results/stage4/ds3-crossed-packet/{label}{ext}' for label in LABELS for ext in ['.json','.npz']]
    paths += ['results/stage4/ds3-crossed-packet/accelerator-validation.json',
              'results/stage4/ds3-crossed-packet/bounds-validation.json']
    summary['hashes']={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}
    (BASE/'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    print(json.dumps({k:summary[k] for k in ['status','Ns','Nt','D','cauchy_ratios','grid_n6_to_n8','same_grid_series_change','total_kernel_points']},indent=2))


if __name__=='__main__':run()
