"""Independent chord polynomials and Fourier validation of the cubic limit."""
from pathlib import Path
import hashlib,json,platform
import sympy as s
import mpmath as mp
from singlet_coefficient import pairings
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-filtered-spectrum'
q,r=s.symbols('q r');lam=s.symbols('lam',positive=True)
def moment(u,v):
    val=0
    for pa in pairings(tuple(range(u+v))):
        cr=sum(a<c<b<d or c<a<d<b for i,(a,b) in enumerate(pa) for c,d in pa[i+1:])
        br=sum(a<u<=b for a,b in pa)
        val+=q**cr*r**br
    return s.expand(val)
def main():
    target=OUT/'spectral-checks.json'
    if target.exists():raise FileExistsError(target)
    entries=[dict(u=u,v=v,polynomial=str(moment(u,v))) for total in [0,2,4] for u in range(total+1) for v in [total-u]]
    assert s.expand(moment(4,0)-(2+q))==0
    assert s.expand(moment(3,1)-r*(2+q))==0
    assert s.expand(moment(2,2)-(1+(1+q)*r*r))==0
    d4=s.expand(sum((-1)**j*s.binomial(4,j)*moment(j,4-j) for j in range(5)))
    assert s.expand(d4-(10+2*q-(16+8*q)*r+6*(1+q)*r*r))==0
    j2d2=s.expand((moment(4,0)+moment(0,4)-2*moment(2,2))/4)
    assert s.expand(j2d2-(1+q)*(1-r*r)/2)==0
    a=-(1-q)**2/(1+q*q);D=(1+q*q)*(1+q+q*q)
    A1=4*(1-q*q)*(1+q*q)*(1+q+q*q)**2
    A2=16*q*q*(1+q)*(1+q*q)*(1+q+q*q)*(2*q*q-q+2)
    VR=(lam*A1+lam*lam*A2)/(lam**3*D**2);VF=64*q*q;W=4*(1/lam+4*q/(1+q))
    norm=VR+VF+W/lam*(a*a*(1+q*q)/2+s.Rational(9,2)*(1-q*q))
    m2=VR*2/lam*(1-q**3)+VF*2/lam*(1-q)+W/lam**2*(a*a*j2d2.subs(r,q*q)+s.Rational(9,4)*d4.subs(r,q*q))
    # First compute the Laurent leading coefficients exactly.
    leadnorm=s.limit(lam*norm.subs(q,s.exp(-4*lam)),lam,0,dir='+')
    leadm2=s.limit(lam*m2.subs(q,s.exp(-4*lam)),lam,0,dir='+')
    assert leadnorm==176 and leadm2==8832
    fourier=[]
    for omega in [0,1,4,10]:
        vals=[];conditions=[]
        for dps in [50,80]:
            mp.mp.dps=dps;w=mp.mpf(omega)
            pieces=[2*mp.quad(lambda t,k=k:mp.cos(w*t)/mp.cosh(2*t)**k,[0,1,2,4,8,16,32,mp.inf]) for k in [4,6]]
            actual=(47*pieces[1]-36*pieces[0])/11
            S2=mp.mpf(2)/3 if w==0 else mp.pi*w*(w*w+16)/(96*mp.sinh(mp.pi*w/4))
            exact=(47*w*w+128)*S2/880
            assert abs(actual/exact-1)<mp.mpf('1e-40')
            vals.append(mp.nstr(actual,46));conditions.append(mp.nstr((47*abs(pieces[1])+36*abs(pieces[0]))/(11*abs(actual)),25))
        assert abs(mp.mpf(vals[0])/mp.mpf(vals[1])-1)<mp.mpf('1e-40')
        fourier.append(dict(omega=omega,values=vals,cancellation_factors=conditions))
    names=['code/check_filtered_spectrum.py','code/singlet_coefficient.py']
    target.write_text(json.dumps(dict(status='passed',chord_moments=entries,
        exact_leading_norm=str(leadnorm),exact_leading_second_moment=str(leadm2),
        normalized_limit_mu2=str(leadm2/leadnorm),fourier=fourier,precisions=[50,80],
        environment=dict(python=platform.python_version(),sympy=s.__version__,mpmath=mp.__version__),
        input_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names}),indent=2)+'\n')
    print('Exact chord moments, Laurent limits, and four Fourier points passed.')
if __name__=='__main__':main()
