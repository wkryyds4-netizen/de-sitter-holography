"""Replay Bose observables from unsimplified E/W/rho and retained complex U."""
import hashlib
import json
from pathlib import Path
import numpy as np
from replay_ds3_crossed_packet import E,K

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'results/stage4/ds3-bose'
LABELS=['pilot-n4-L8','fixed-n4-L12','refined-n6-L12']


def run():
    results=[]
    for label in LABELS:
        data=json.loads((BASE/(label+'.json')).read_text());assert data['status']=='completed'
        path=ROOT/data['raw_grid'];assert hashlib.sha256(path.read_bytes()).hexdigest()==data['raw_grid_sha256']
        with np.load(path) as z:
            x,w,cs,cw,U,cond=[z[k] for k in ['x','w','c','cw','U','condition']]
        a,b,d=np.meshgrid(x,x,x,indexing='ij');a,b,d=[v.ravel() for v in (a,b,d)]
        ww=(w[:,None,None]*w[None,:,None]*w[None,None,:]).ravel()
        ts=data['times'];n=0.;Y=np.zeros(len(ts),complex);E1=Y.copy();E2=Y.copy();Z=Y.copy();dx=Y.copy();dd=np.zeros(len(ts))
        worst=[];weighted_hermiticity=0.;min_surface=float('inf')
        for i,(c,wc) in enumerate(zip(cs,cw)):
            Wb=2*b*b/np.pi*K(a,b,.6)*K(b,c,.6);Wd=2*d*d/np.pi*K(a,d,.6)*K(d,c,.6)
            amp=E(a,b,.6)*E(b,c,.6)*E(a,d,.6)*E(d,c,.6)/np.sqrt(Wb*Wd)/(2*c*c/np.pi)
            amp*=np.exp(np.pi*(np.sqrt(1+a*a)-np.sqrt(1+c*c)))
            wt=wc*ww*amp
            flip=U[i].reshape((len(x),)*3).swapaxes(1,2).conj().ravel()
            weighted_hermiticity+=np.sum(wt*abs(U[i]-flip))
            proximity=abs(a+c-b-d);min_surface=min(min_surface,float(proximity.min()))
            sensitivity=wt*abs(U[i])*cond[i]*np.finfo(float).eps
            for j in np.argpartition(sensitivity,-3)[-3:]:
                worst.append(dict(a=float(a[j]),b=float(b[j]),c=float(c),d=float(d[j]),
                                  weighted_proxy=float(sensitivity[j]),weight=float(wt[j]),
                                  condition=float(cond[i,j]),u=[float(U[i,j].real),float(U[i,j].imag)]))
            hA,hB,hD,hC=np.sqrt(1+a*a),np.sqrt(1+b*b),np.sqrt(1+d*d),np.sqrt(1+c*c)
            for k,t in enumerate(ts):
                alpha=np.exp(-1j*t*(hA-hB));beta=np.exp(-1j*t*(hB-hC))
                alphaD=np.exp(-1j*t*(hA-hD));betaD=np.exp(-1j*t*(hD-hC))
                E1[k]+=np.sum(wt*U[i]*alpha.conj()*alphaD)
                E2[k]+=np.sum(wt*U[i]*beta.conj()*betaD)
                Z[k]+=np.sum(wt*U[i]*alpha.conj()*betaD)
                dx[k]+=np.sum(wt*U[i]*(alpha-beta).conj()*(alphaD-betaD))
            aa=x[:,None];bb=x[None,:]
            W=2*bb*bb/np.pi*K(aa,bb,.6)*K(bb,c,.6)
            w2=wc*w[:,None]*w[None,:]*E(aa,bb,.6)**2*E(bb,c,.6)**2/W/(2*c*c/np.pi)
            w2*=np.exp(np.pi*(np.sqrt(1+aa*aa)-hC))
            n+=w2.sum()
            for k,t in enumerate(ts):
                al=np.exp(-1j*t*(np.sqrt(1+aa*aa)-np.sqrt(1+bb*bb)))
                be=np.exp(-1j*t*(np.sqrt(1+bb*bb)-hC))
                Y[k]+=np.sum(w2*al.conj()*be);dd[k]+=np.sum(w2*abs(al-be)**2)
        residual=[]
        for k,r in enumerate(data['rows']):
            for name,calc in [('E1',E1[k]),('E2',E2[k]),('Y',Y[k]),('Z',Z[k]),('D_exchange',dx[k]),('D',dd[k]+dx[k])]:
                residual.append(abs(calc-complex(*r[name])))
            residual.extend([abs(n-r['N_label']),abs(dd[k]-r['D_direct'])])
            assert r['N1'][0]>0 and r['N2'][0]>0 and r['D'][0]>=-1e-17
            assert r['cauchy_ratio']<=1+1e-12
            assert abs(r['expanded_D']-r['D'][0])<1e-17
        assert max(residual)<3e-17
        assert data['rows'][0]['D']==[0.,0.]
        results.append(dict(label=label,max_replay_difference=float(max(residual)),
                            weighted_hermiticity=float(weighted_hermiticity),min_a_plus_c_minus_b_minus_d=min_surface,
                            worst_weighted_nodes=sorted(worst,key=lambda r:r['weighted_proxy'],reverse=True)[:6],
                            points=data['points']))
    data=[json.loads((BASE/(label+'.json')).read_text()) for label in LABELS]
    differences=[]
    for i,j,meaning in [(0,1,'L and initial-node change; not isolated tail'),(1,2,'only mass-node order 4 to 6; fixed L12 nc12 series128 K12')]:
        a=data[i];b=data[j]
        differences.append(dict(comparison=meaning,D=[abs(complex(*x['D'])-complex(*y['D'])) for x,y in zip(a['rows'],b['rows'])],
                                overlap=[abs(complex(*x['overlap'])-complex(*y['overlap'])) for x,y in zip(a['rows'],b['rows'])]))
    paths=['code/replay_ds3_bose.py','code/replay_ds3_crossed_packet.py']+[f'results/stage4/ds3-bose/{l}{ext}' for l in LABELS for ext in ['.json','.npz']]
    summary=dict(status='passed',results=results,differences=differences,total_points=sum(r['points'] for r in results),
                 final_label=LABELS[-1],error_scope='Empirical replay/refinement and identities; not total error certification',
                 hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths})
    (BASE/'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    print(json.dumps({k:summary[k] for k in ['status','differences','total_points']},indent=2))


if __name__=='__main__':run()
