"""Analytic norm bounds for compact initial probability and genuine infinite tails."""
import mpmath as mp


def bounds(L,s='.2',c0='.4',c1='.8',p='.6',q='.9'):
    L,s,c0,c1,p,q=map(mp.mpf,[L,s,c0,c1,p,q]);pi=mp.pi
    assert 0<s<=c0<c1<L
    delta=lambda x:1/(mp.sqrt(1+x*x)+x)
    Cs=mp.exp(pi*delta(s));CL=mp.exp(pi*delta(L))
    gmin=c0*mp.exp(pi*delta(c1))/(1+mp.exp(-2*pi*c0))
    decay=mp.exp(-pi*(L-c1))
    tail_a=CL*decay/(16*pi**3*gmin)*(L*L+L*(3/pi-c0)+3/pi**2-c0/pi)
    C0=(Cs-1)/(2*delta(s))
    tail_b=decay/(8*pi**4*gmin)*(L+C0+2/pi)
    tail=tail_a+tail_b
    Ap=p/(4*pi*mp.sinh(pi*p));Aq=q/(4*pi*mp.sinh(pi*q))
    moment=1/(2*pi**3)
    norm=((c1+C0)*Ap*Aq+(Ap+Aq)*moment)/gmin
    return dict(L=L,gmin=gmin,Cs=Cs,C0=C0,tail_a=tail_a,tail_b=tail_b,
                tail_norm_bound=tail,norm_squared_bound=norm,
                crossed_absolute_bound=2*mp.sqrt(tail*norm))
