"""Gamma-phase certificates for the physical collision peak (not a mass integral)."""
import sympy as sp


def gamma_terms(b,d,a,q,c,p,kind):
    hs=[(1+sp.I*x)/2 for x in [a,q,c,p]];h1,h2,h3,h4=hs;h=(1+sp.I*b)/2;hp=(1+sp.I*d)/2
    terms=[]
    def add(top,bottom,w=1):
        terms.extend((sp.expand(x),w) for x in top);terms.extend((sp.expand(x),-w) for x in bottom)
    def shadow(H,hi,hj):
        dd=hi-hj;add([2*H-1,1-H-dd,1-H+dd],[2-2*H,H-dd,H+dd])
    def pref():add([h+h1-h2,h-h1+h2,1-h+h4-h3,1-h+h3-h4],[2*h,2-2*h])
    def reflect(H,Q):add([H+h1+h2-1,H+h3+h4-1,Q+1-h1-h4,Q+1-h2-h3],[H+1-h1-h2,H+1-h3-h4,Q+h1+h4-1,Q+h2+h3-1])
    def period(nus,w):
        S=sum(nus);add([(1-sp.I*S)/2]+[(1-sp.I*(S-2*x))/2 for x in nus],[1-sp.I*x for x in nus],w)
    shadow(1-hp,h1,h4);pref()
    if kind=='conjugate_coefficient':
        shadow(hp,h3,h2);pref();reflect(h,1-hp);reflect(1-h,1-hp);weight=2
    elif kind=='peak_phase':
        z=h-hp+h1+h3-1
        A=[h+h1+h2-1,h+h3+h4-1,h+h3-h4,h+h1-h2]
        add([z,2*h,A[2]-z,A[3]-z],[A[2],A[3],2*h-z]);weight=1
    else:raise ValueError(kind)
    period([a,q,b],weight);period([c,p,-b],weight);period([c,q,d],-weight);period([a,p,-d],-weight)
    counts={}
    for arg,w in terms:counts[arg]=counts.get(arg,0)+w
    return {arg:w for arg,w in counts.items() if w}


def symbolic_certificates():
    b,d,a,q,c,p=sp.symbols('b d a q c p',real=True)
    counts=gamma_terms(b,d,a,q,c,p,'conjugate_coefficient')
    pairs=all(w==counts.get(sp.expand(sp.conjugate(arg)),0) for arg,w in counts.items())
    M,u,v=sp.symbols('M u v',real=True);L=sp.Symbol('L',real=True)
    peak=gamma_terms(M+u,M+v,M+u+v,q,M,p,'peak_phase');small={};large={}
    for arg,w in peak.items():
        slope=sp.expand(sp.im(arg)).coeff(M)
        (small if slope==0 else large)[arg]=w
    small_cancel=all((sp.im(arg)==0 and sp.re(arg)>0) or w==small.get(sp.expand(sp.conjugate(arg)),0) for arg,w in small.items())
    phase=0
    for arg,w in large.items():
        slope=sp.expand(sp.im(arg)).coeff(M);aa=sp.expand(arg-sp.I*slope*M)
        logbase=L+sp.log(abs(slope))+sp.I*sp.sign(slope)*sp.pi/2
        expression=(sp.I*slope*M+aa-sp.Rational(1,2))*logbase-sp.I*slope*M
        for n in [1,2]:expression+=(-1)**(n+1)*sp.bernoulli(n+1,aa)/(n*(n+1)*(sp.I*slope*M)**n)
        phase+=w*sp.im(sp.expand(expression))
    phase=sp.factor(phase)
    expected=(-p*p*v+q*q*u+u*u*v-u*v*v+u-v)/(4*M*M)
    return dict(coefficient_conjugate_pairs=pairs,small_phase_cancels=bool(small_cancel),phase_expansion_correct=sp.simplify(phase-expected)==0,
                phase_through_inverse_square=str(phase),coefficient_terms=[dict(argument=str(k),weight=w) for k,w in counts.items()],
                peak_large_terms=[dict(argument=str(k),weight=w) for k,w in large.items()],
                peak_small_terms=[dict(argument=str(k),weight=w) for k,w in small.items()])
