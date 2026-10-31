"""LPRS Omega at a removable Gamma-pair collision, with mp zeta-tail sums.

No claim of uniform asymptotic remainder certification. Exact collision
uses a convergent digamma-weighted series, not displacement of the input.
"""
import mpmath as mp
from global_racah import omega_hyper, omega_quad, shadow_coefficient
from ds3_crossing_dictionary import period


def collision_data(h,hs):
    h1,h2,h3,h4=hs
    hp0=h+h1+h3-1
    top=[h+h1+h2-1,h+h3+h4-1,h+h3-h4,h+h1-h2]
    bottom=[2*hp0,2*h,mp.mpf(1)]
    return hp0,top,bottom


def omega_collision(h,hs,N=128,K=16):
    hp0,top,bottom=collision_data(h,hs)
    term=mp.mpc(1);values=[]
    for n in range(N):
        bracket=mp.digamma(2*h+n)+mp.digamma(2*hp0+n)+2*mp.digamma(1+n)-sum(mp.digamma(a+n) for a in top)
        values.append(term*bracket)
        term*=mp.fprod(a+n for a in top)/(mp.fprod(b+n for b in bottom)*(n+1))
    # Gamma-ratio expansion of t_n=C*n^-2*sum c_k/n^k.
    logs=[mp.mpc(0)];cs=[mp.mpc(1)]
    for k in range(1,K+1):
        diff=sum(mp.bernpoly(k+1,a) for a in top)-sum(mp.bernpoly(k+1,b) for b in bottom)-mp.bernpoly(k+1,1)
        logs.append((-1)**(k+1)*diff/(k*(k+1)))
        cs.append(sum(j*logs[j]*cs[k-j] for j in range(1,k+1))/k)
    ds=[mp.mpc(0)]
    for k in range(1,K+1):
        ds.append((-1)**(k+1)/k*(mp.bernpoly(k,2*h)+mp.bernpoly(k,2*hp0)+2*mp.bernpoly(k,1)-sum(mp.bernpoly(k,a) for a in top)))
    coeff=[sum(cs[j]*ds[k-j] for j in range(k)) for k in range(1,K+1)]
    C=mp.fprod(mp.gamma(b) for b in bottom)/mp.fprod(mp.gamma(a) for a in top)
    tails=[C*v*mp.zeta(k+2,N) for k,v in enumerate(coeff,1)]
    total=mp.fsum(values)+mp.fsum(tails)
    return total,dict(N=N,K=K,terms_absolute=mp.fsum(abs(v) for v in values),
                      tail=mp.fsum(tails),last_tail_term=tails[-1],leading_bracket=ds[1])


def omega_collision_derivative(h,hs):
    """Independent parameter derivative of the two regular 4F3 prefactors."""
    h1,h2,h3,h4=hs;p=h2+h3;R=h1-h2;S=h3-h4
    hp0,top,bottom=collision_data(h,hs)
    def pieces(z):
        hp=hp0-z
        A=mp.gamma(2*h)*mp.gamma(hp-p+1)*mp.gamma(hp-R+S-p+1)/(mp.gamma(R+h)*mp.gamma(S+h)*mp.gamma(hp-R+h-p+1))
        A*=mp.hyper([hp+h2-h3,hp-h1+h4,hp-R+S-p+1,hp-p+1],
                    [2*hp,hp-R+h-p+1,hp-R-h-p+2],1)
        B=mp.gamma(2*hp)*mp.gamma(h1-h3+h+p-1)*mp.gamma(h4-h2+h+p-1)/(mp.gamma(hp+h2-h3)*mp.gamma(hp-h1+h4)*mp.gamma(hp+R+h+p-1))
        B*=mp.hyper([h1-h3+h+p-1,h4-h2+h+p-1,S+h,R+h],
                    [hp+R+h+p-1,2*h,-hp+R+h+p],1)
        return A,B
    f0=mp.hyper(top,bottom,1)
    return mp.diff(lambda z:pieces(z)[0]-pieces(z)[1],0)-2*mp.euler*f0


def physical_resolved(b,d,external,N=128,K=16,method='series'):
    """Exact collision only; generic points use high-precision LPRS hyper sums.

    Zero detection tolerates context roundoff, not an O(1e-5) mass window.
    """
    b,d=mp.mpf(b),mp.mpf(d);external=list(map(mp.mpf,external))
    hs=[(1+1j*x)/2 for x in external];h,hp=(1+1j*b)/2,(1+1j*d)/2
    h1,h2,h3,h4=hs
    pref=mp.gamma(h+h1-h2)*mp.gamma(h-h1+h2)*mp.gamma(1-h+h4-h3)*mp.gamma(1-h+h3-h4)/(4*mp.gamma(2*h)*mp.gamma(2-2*h))
    terms=[];oms=[]
    for q,co in [(hp,shadow_coefficient(1-hp,h1,h4)),(1-hp,shadow_coefficient(hp,h3,h2))]:
        vs=[]
        for H in [h,1-h]:
            z=H-q+h1+h3-1
            hit=abs(z)<64*mp.eps*max(1,abs(H),abs(q),*[abs(x) for x in hs])
            if method=='euler':val=omega_quad(H,q,hs);meta={}
            elif hit:
                if method=='derivative':val=omega_collision_derivative(H,hs);meta={}
                else:val,meta=omega_collision(H,hs,N,K)
            else:val,cond=omega_hyper(H,q,hs);meta=dict(condition=cond)
            oms.append(dict(collision=hit,z=z,value=val,**meta));vs.append(val)
        terms.append(co*pref*mp.sin(mp.pi*(q-h1-h4))*mp.sin(mp.pi*(q-h2-h3))*vs[0]*vs[1])
    a,q,c,p=external;ns=period(a,q,b)*period(c,p,-b);nt=period(c,q,d)*period(a,p,-d)
    factor=b*d/mp.pi**4*(ns/abs(ns))/(nt/abs(nt))
    value=factor*mp.fsum(terms)
    return value,dict(omega=oms,shadow_condition=mp.fsum(abs(v) for v in terms)/abs(mp.fsum(terms)))
