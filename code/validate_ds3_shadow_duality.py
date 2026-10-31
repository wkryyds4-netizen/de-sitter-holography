"""Symbolic dictionary, saved original factors, generic hyper and independent Euler checks."""
import hashlib,json,platform,time
from pathlib import Path
import mpmath as mp
import sympy as sp
from ds3_shadow_duality import reflection,wilson_data
from global_racah import omega_hyper,omega_quad
from run_ds3_collision import clean


def symbolic():
    masses=sp.symbols('b d a q c p',real=True)
    b,d,a,q,c,p=masses;hs=[(1+sp.I*v)/2 for v in [a,q,c,p]];H=(1+sp.I*b)/2;hp=(1+sp.I*d)/2
    def data(H,hp):
        h1,h2,h3,h4=hs;hp0=H+h1+h3-1;z=hp0-hp
        A=[H+h1+h2-1,H+h3+h4-1,H+h3-h4,H+h1-h2]
        return [sp.expand(v-z/2) for v in A],sp.expand(z/(2*sp.I)),sp.expand((hp0-H-z/2)/sp.I)
    T,l,x=data(H,hp);Ts,ls,xs=data(1-H,1-hp)
    checks={'lambda_to_x':sp.simplify(l-xs)==0,'x_to_lambda':sp.simplify(x-ls)==0,
            'dual_to_conjugate':all(sp.simplify(sum(T)/2-t-sp.conjugate(s))==0 for t,s in zip(T,Ts))}
    return dict(checks=checks,T=list(map(str,T)),lam=str(l),x=str(x),reflected_T=list(map(str,Ts)))


def main():
    root=Path(__file__).resolve().parents[1];base=root/'results/stage4/ds3-shadow-duality';base.mkdir(exist_ok=True)
    target=base/'validation.json';raw=base/'validation.jsonl'
    if raw.exists() or target.exists():raise RuntimeError('Preserve prior checks')
    rows=[]
    def save(name,got,ref,tol,extra=None):
        error=abs(got-ref)/max(abs(ref),mp.mpf('1e-1000'))
        r=clean(dict(name=name,value=got,reference=ref,relative_error=error,tolerance=tol,passed=error<mp.mpf(tol),extra=extra))
        rows.append(r)
        with raw.open('a',encoding='utf-8') as f:f.write(json.dumps(r)+'\n');f.flush()
        print(name,mp.nstr(error,5),r['passed'],flush=True)
    mp.mp.dps=80
    old=json.loads((root/'results/stage4/ds3-collision/study.json').read_text())
    for row in old['rows']:
        M=mp.mpf(row['M']);r=mp.mpf(row['offset']);a=M+mp.mpf('.1');b=M+mp.mpf('.3');d=M-mp.mpf('.2')+r
        hs=[(1+1j*v)/2 for v in [a,mp.mpf('.6'),M,mp.mpf('.6')]];H=(1+1j*b)/2;hp=(1+1j*d)/2
        O=[mp.mpc(*v['value']) for v in row['details'][-1]['omega']]
        for h,i,j in [(H,2,1),(1-H,3,0)]:
            pred=reflection(h,1-hp,hs)*mp.conj(O[j])
            save(f'saved_M{M}_r{r}_factor{i}',pred,O[i],'1e-24' if r==0 else '1e-60',{'prototype_bias_retained':r==0})
    for index,labels in enumerate([['.7','.9','1.1','.4','.8','.6'],['1.3','.8','.9','.7','1.6','.3']]):
        mp.mp.dps=50;b,d,a,q,c,p=map(mp.mpf,labels);H=(1+1j*b)/2;hp=(1+1j*d)/2;hs=[(1+1j*v)/2 for v in [a,q,c,p]]
        for h in [H,1-H]:
            direct=omega_hyper(h,1-hp,hs)[0]
            reflected=reflection(h,1-hp,hs)*mp.conj(omega_hyper(1-h,hp,hs)[0])
            save(f'generic_hyper_{index}_{mp.nstr(h,5)}',reflected,direct,'1e-40',{'labels':labels})
            K=wilson_data(h,1-hp,hs)[-1];Ks=wilson_data(1-h,hp,hs)[-1]
            save(f'gamma_ratio_{index}_{mp.nstr(h,5)}',reflection(h,1-hp,hs),K/mp.conj(Ks),'1e-45')
    mp.mp.dps=35;b,d,a,q,c,p=map(mp.mpf,['.7','.9','1.1','.4','.8','.6']);H=(1+1j*b)/2;hp=(1+1j*d)/2;hs=[(1+1j*v)/2 for v in [a,q,c,p]]
    direct=omega_quad(H,1-hp,hs);other=omega_quad(1-H,hp,hs)
    save('independent_Euler_pair',reflection(H,1-hp,hs)*mp.conj(other),direct,'1e-14',{'dps':35,'method':'two direct Euler product integrals'})
    sym=symbolic()
    paths=['code/ds3_shadow_duality.py','code/validate_ds3_shadow_duality.py','code/global_racah.py','results/stage4/ds3-collision/study.json']
    result=dict(status='passed' if all(r['passed'] for r in rows) and all(sym['checks'].values()) else 'failed',checks=rows,symbolic=sym,
                versions=dict(python=platform.python_version(),mpmath=mp.__version__,sympy=sp.__version__),
                hashes={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in paths})
    target.write_text(json.dumps(result,indent=2),encoding='utf-8');print(result['status'],len(rows),flush=True)


if __name__=='__main__':main()
