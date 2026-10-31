"""Scalar d=2 conformal 6j, LPRS 1808.00612v3 (3.17),(3.35)-(3.38).

hbar=h throughout (scalar principal series), NOT hbar=conjugate(h).
This is a matrix-element evaluator, not a full spin-summed transform.
"""
import mpmath as mp


def omega_hyper(h, hp, hs):
    h1,h2,h3,h4=hs; p=h2+h3
    h12,h34,h23,h14,h13,h42=h1-h2,h3-h4,h2-h3,h1-h4,h1-h3,h4-h2
    a=(mp.gamma(2*h)*mp.gamma(hp-p+1)*mp.gamma(hp-h12+h34-p+1)
       *mp.gamma(-hp+h12+h+p-1)/(mp.gamma(h12+h)*mp.gamma(h34+h)
                               *mp.gamma(hp-h12+h-p+1)))
    a*=mp.hyper([hp+h23,hp-h14,hp-h12+h34-p+1,hp-p+1],
                [2*hp,hp-h12+h-p+1,hp-h12-h-p+2],1)
    b=(mp.gamma(2*hp)*mp.gamma(hp-h12-h-p+1)*mp.gamma(h13+h+p-1)
       *mp.gamma(h42+h+p-1)/(mp.gamma(hp+h23)*mp.gamma(hp-h14)
                            *mp.gamma(hp+h12+h+p-1)))
    b*=mp.hyper([h13+h+p-1,h42+h+p-1,h34+h,h12+h],
                [hp+h12+h+p-1,2*h,-hp+h12+h+p],1)
    value=a+b
    return value, (abs(a)+abs(b))/abs(value)


def omega_quad(h,hp,hs):
    """Direct, convergent Euler product integral (3.37) for our samples."""
    h1,h2,h3,h4=hs; p=h2+h3
    def f(x):
        return (x**(h+p+h1-h3-2)*(1-x)**(hp-p)
                *mp.hyp2f1(h+h1-h2,h-h3+h4,2*h,x)
                *mp.hyp2f1(hp-h3+h2,hp+h1-h4,2*hp,1-x))
    return mp.quad(f,[0,mp.mpf('.01'),mp.mpf('.5'),mp.mpf('.99'),1])


def shadow_coefficient(h,hi,hj):
    d=hi-hj
    return (mp.pi*mp.gamma(2*h-1)/mp.gamma(2-2*h)
            *mp.gamma(1-h-d)*mp.gamma(1-h+d)/(mp.gamma(h-d)*mp.gamma(h+d)))


def scalar_6j(s,t,external,method='hyper'):
    """s,t are principal-series nu labels; external=(nu1,nu2,nu3,nu4)."""
    hs=[(1+1j*x)/2 for x in external]
    h,hp=(1+1j*s)/2,(1+1j*t)/2
    h1,h2,h3,h4=hs
    pref=(mp.gamma(h+h1-h2)*mp.gamma(h-h1+h2)
          *mp.gamma(1-h+h4-h3)*mp.gamma(1-h+h3-h4)
          /(4*mp.gamma(2*h)*mp.gamma(2-2*h)))
    om=[]; condition=[]
    def B(q):
        vals=[]
        for x in [h,1-h]:
            if method=='quad':
                val=omega_quad(x,q,hs); cond=None
            else:
                val,cond=omega_hyper(x,q,hs)
            vals.append(val);om.append(val);condition.append(cond)
        return pref*mp.sin(mp.pi*(q-h1-h4))*mp.sin(mp.pi*(q-h2-h3))*vals[0]*vals[1]
    terms=[shadow_coefficient(1-hp,h1,h4)*B(hp),
           shadow_coefficient(hp,h3,h2)*B(1-hp)]
    val=sum(terms)
    return dict(value=val,omega=om,omega_cancellation=condition,
                shadow_terms=terms,shadow_cancellation=sum(abs(x) for x in terms)/abs(val))


def spectral_rectangle(values):
    """Row-major 2x2 invariant under nonzero row/column rescalings."""
    return values[0]*values[3]/(values[1]*values[2])


def shadow_fourier(nu,p):
    """Unnormalized d=2 scalar shadow multiplier; FT uses exp(-ip.x)."""
    return mp.pi*2**(2j*nu)*mp.gamma(1j*nu)/mp.gamma(1-1j*nu)*p**(-2j*nu)
