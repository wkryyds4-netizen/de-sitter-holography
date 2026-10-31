"""Near-collision Omega: direct sums with analytic absolute series-tail bounds.

Bounds concern exact arithmetic. mp precision refinements separately
check rounding; this is not interval arithmetic or a full-6j bound.
"""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from ds3_collision_asymptotic import leading_bound
from run_ds3_collision import clean


def hyper_sum(top,bottom,N,m):
    t=mp.mpc(1);values=[];product=mp.mpf(1)
    for n in range(N):
        values.append(t)
        t*=mp.fprod(n+a for a in top)/(mp.fprod(n+b for b in bottom)*(n+1))
        product*=mp.mpf(n+1)**2/(mp.mpf(n+1)**2+m*m)
    return mp.fsum(values),product/N


def controlled(M,u,v,p,q,r='0',N=128):
    M,u,v,p,q,r=map(mp.mpf,[M,u,v,p,q,r])
    b=M+u;d0=M+v;d=d0+r;H=(1-1j*b)/2;hp0=(1+1j*d0)/2;z=-1j*r/2
    top=[(1+1j*(v+q))/2,(1+1j*(p-u))/2,(1-1j*(u+p))/2,(1+1j*(v-q))/2]
    # Check the theorem's hypotheses, including all shifted numerator parameters.
    if min(b,d,d0)>0 and all(abs(mp.im(a))<=mp.mpf('.5')+32*mp.eps for a in top) and all(abs(mp.im(a-z))<=mp.mpf('.5')+32*mp.eps for a in top):pass
    else:raise ValueError('Outside proved numerator/mass window')
    if r==0:
        lead,rem,detail=leading_bound(M,u,v,p,q)
        t=mp.mpc(1);values=[];product=mp.mpf(1);m=min(b,d0)
        bottom=[2*hp0,2*H,mp.mpf(1)]
        for n in range(N):
            bracket=mp.digamma(2*H+n)+mp.digamma(2*hp0+n)+2*mp.digamma(1+n)-sum(mp.digamma(a+n) for a in top)
            values.append(t*bracket)
            t*=mp.fprod(n+a for a in top)/(mp.fprod(n+x for x in bottom)*(n+1))
            product*=mp.mpf(n+1)**2/(mp.mpf(n+1)**2+m*m)
        value=mp.fsum(values);tail=detail['bracket_bound']*product/N
        return value,dict(leading=lead,leading_bound=rem,series_tail_bound=tail,N=N,collision=True,
                          summation_condition=mp.fsum(abs(x) for x in values)/abs(value))
    hp=hp0-z
    prefA=mp.gamma(2*H)*mp.gamma(top[3]-z)*mp.gamma(top[2]-z)/(mp.gamma(top[3])*mp.gamma(top[2])*mp.gamma(2*H-z))
    prefB=mp.gamma(2*hp)*mp.gamma(top[0])*mp.gamma(top[1])/(mp.gamma(top[0]-z)*mp.gamma(top[1]-z)*mp.gamma(2*hp0-z))
    ca=mp.gamma(z)*prefA;cb=mp.gamma(-z)*prefB
    ma=min(abs(d),abs(b-r/2));mb=min(abs(d0+r/2),abs(b))
    if min(ma,mb)<=0:raise ValueError('Degenerate mass denominator')
    fa,ta=hyper_sum([a-z for a in top],[2*hp,2*H-z,1-z],N,ma)
    fb,tb=hyper_sum(top,[2*hp0-z,2*H,1+z],N,mb)
    va,vb=ca*fa,cb*fb;value=va+vb
    remainder=(mp.zeta(2)-1)*(abs(ca)/(1+ma*ma)+abs(cb)/(1+mb*mb))
    return value,dict(leading=ca+cb,leading_bound=remainder,series_tail_bound=abs(ca)*ta+abs(cb)*tb,
                      N=N,collision=False,summation_condition=(abs(va)+abs(vb))/abs(value),mA=ma,mB=mb)


def main():
    root=Path(__file__).resolve().parents[1];dest=root/'results/stage4/ds3-collision/controlled.json'
    if dest.exists():raise RuntimeError('Preserve prior run')
    rows=[]
    for M in [16,32,64,256,1024]:
        for lam in [0,1]:
            vals=[];records=[]
            for dps in [70,100]:
                with mp.workdps(dps):
                    r=mp.mpf(lam)/mp.log(M) if lam else mp.mpf(0)
                    val,detail=controlled(M,'.3','-.2','.6','.6',r,128)
                    vals.append(val);records.append(clean(detail))
            with mp.workdps(100):
                target=mp.sin(mp.mpf(lam)/2)/(mp.mpf(lam)/2) if lam else mp.mpf(1)
                scaled=vals[-1]/(2*mp.log(M))
                rows.append(dict(M=M,lambda_value=lam,r=str(r),value=clean(vals[-1]),
                                 scaled=clean(scaled),sinc_target=str(target),scaled_error=str(abs(scaled-target)),
                                 precision_difference=str(abs(vals[0]-vals[1])),details=records))
    paths=['code/ds3_controlled_omega.py','code/ds3_collision_asymptotic.py','code/ds3_racah_merged.py','code/run_ds3_collision.py']
    data=dict(status='completed',precisions=[70,100],rows=rows,
              scope='Exact-arithmetic analytic series-tail bounds; no interval rounding or full 6j certification',
              hashes={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in paths})
    dest.write_text(json.dumps(data,indent=2),encoding='utf-8')
    print('completed',len(rows))


if __name__=='__main__':main()
