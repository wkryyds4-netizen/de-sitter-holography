"""Explicit collision-Omega bound on |u|+p<=1, |v|+q<=1.

This is an Omega-factor bound, not the error of the complete 6j or Bose integral.
"""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from ds3_racah_merged import omega_merged
from run_ds3_collision import clean


def leading_bound(M,u,v,p,q):
    M,u,v,p,q=map(mp.mpf,[M,u,v,p,q])
    if abs(u)+p>1 or abs(v)+q>1 or min(M+u,M+v)<=0:
        raise ValueError('Outside proved transfer window')
    b=M+u;d=M+v
    top=[(1+1j*(v+q))/2,(1+1j*(p-u))/2,
         (1-1j*(u+p))/2,(1+1j*(v-q))/2]
    leading=mp.digamma(1-1j*b)+mp.digamma(1+1j*d)+2*mp.digamma(1)-sum(mp.digamma(a) for a in top)
    bracket_bound=2+mp.asinh(abs(b))+mp.asinh(abs(d))+4*mp.log(2)*sum(abs(a-1) for a in top)
    remainder=(mp.zeta(2)-1)*bracket_bound/(1+min(b,d)**2)
    return leading,remainder,dict(top=top,bracket_bound=bracket_bound,b=b,d=d)


def main():
    root=Path(__file__).resolve().parents[1];dest=root/'results/stage4/ds3-collision/asymptotic.json'
    if dest.exists():raise RuntimeError('Preserve prior run')
    mp.mp.dps=70;rows=[]
    for u,v,p,q in [('.3','-.2','.6','.6'),('.4','-.3','.6','.7')]:
        for M in [2,4,8,16,32,64]:
            mm,uu,vv,pp,qq=map(mp.mpf,[M,u,v,p,q]);a=mm+uu+vv;b=mm+uu;c=mm;d=mm+vv
            H=(1-1j*b)/2;hp=(1+1j*d)/2;hs=[(1+1j*x)/2 for x in [a,qq,c,pp]]
            val,detail=omega_merged(H,hp,hs,max(256,16*M),32)
            lead,bound,terms=leading_bound(M,u,v,p,q)
            err=abs(val-lead)
            rows.append(dict(M=M,u=u,v=v,p=p,q=q,value=val,leading=lead,absolute_error=err,
                             bound=bound,error_to_bound=err/bound,passed=err<bound,details=terms))
    paths=['code/ds3_collision_asymptotic.py','code/ds3_racah_merged.py','code/ds3_racah_collision.py']
    data=dict(status='passed' if all(r['passed'] for r in rows) else 'failed',dps=70,rows=rows,
              scope='The theorem is analytic; evaluations are empirical checks inside its specified window. Not a full 6j limit.',
              hashes={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in paths})
    dest.write_text(json.dumps(clean(data),indent=2),encoding='utf-8')
    print(data['status'],len(rows))


if __name__=='__main__':main()
