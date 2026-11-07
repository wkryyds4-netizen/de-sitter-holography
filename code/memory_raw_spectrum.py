"""Restore raw quartic descendants and the conserved-subtracted memory spectrum."""
from pathlib import Path
import json,hashlib,platform
from functools import lru_cache
import sympy as s

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-memory-raw-spectrum'
q,lam,C,x,y,r=s.symbols('q lam C x y r')


def coefficients():
    old=json.loads((ROOT/'results/susskind-memory-descendants/symbolic-filters.json').read_text())
    parse=lambda z:s.sympify(z.replace('lambda','lam'),locals={'q':q,'lam':lam})
    raw=list(map(parse,old['raw_Q4_coefficients']))
    alpha,beta,gamma,delta=[-raw[i] for i in [0,1,2,4]]
    table={k:parse(v) for k,v in old['raw_Q4_combination'].items()}
    a=-(1-q)**2/(1+q*q)
    gJ=alpha+beta+gamma-1;bD=(alpha-beta)/4;bR=gamma-2
    bJ=gJ+a*bR;bF=bR+table['JhF_coefficient']
    v=2*lam*(1-q)+16*lam**2*q/(1+q);d=2*lam*(1+q);rho=v/(v+d)
    kappa=table['O_coefficient']+8*lam*q*q*rho
    assert s.factor(bF-2*(q-1))==0
    assert s.factor(bD-q*q/(4*sum(q**j for j in range(5))))==0
    assert s.factor(table['O_coefficient']-(1-q)*(1+q)**2/sum(q**j for j in range(5)))==0
    return dict(alpha=alpha,beta=beta,gamma=gamma,delta=delta,a=a,gJ=gJ,
                bD=bD,bR=bR,bJ=bJ,bF=bF,kappa=kappa,rho=rho,**table)


@lru_cache(None)
def monic_power(n):
    v=[s.Integer(1)]
    for k in range(n):
        w=[s.Integer(0)]*(len(v)+1)
        for j,c in enumerate(v):
            w[j+1]+=c
            if j:w[j-1]+=c*sum(q**a for a in range(j))
        v=list(map(s.expand,w))
    return v


def bivariate_moment(a,b):
    va,vb=monic_power(a),monic_power(b);fact=s.Integer(1);out=0
    for n in range(min(a,b)+1):
        if n:fact*=sum(q**j for j in range(n))
        out+=va[n]*vb[n]*fact*r**n
    return s.expand(out)


def average(poly):
    return s.expand(sum(co*bivariate_moment(i,j) for (i,j),co in s.Poly(s.expand(poly),x,y).terms()))


def main():
    OUT.mkdir(exist_ok=True);target=OUT/'symbolic-spectrum.json'
    if target.exists():raise FileExistsError(target)
    co=coefficients();e=(x+y)/2;om=x-y
    # Exact free-algebra raw identity expressed in its five-word coordinates.
    p,m=s.symbols('p m')
    Z=s.Matrix([0,0,0,0,1]);E=s.Matrix([0,0,s.Rational(1,2),0,2*p-m])
    M=s.Matrix([0,0,0,1,2*p-m]);S=s.Matrix([0,s.Rational(1,2),-s.Rational(1,2),-1,m-p])
    CC=s.Matrix([1,-1,0,1,0]);core=CC-2*S+(E+M)/2
    raw=CC-co['alpha']*E-co['beta']*M-co['gamma']*S-co['delta']*Z
    J2=(E+M)/2;L2=2*(E-M);JT3=S-J2
    diff=core-raw-co['gJ']*J2-co['bD']*L2-co['bR']*JT3-co['delta']*Z
    assert all(s.factor(z)==0 for z in diff)
    powers={'e2':e**2,'omega2':om**2,'e4':e**4,'e2omega2':e**2*om**2,'omega4':om**4,
            'e4omega2':e**4*om**2,'e2omega4':e**2*om**4,'omega6':om**6}
    moments={k:s.factor(average(v)) for k,v in powers.items()}
    assert s.factor(moments['e2']-(1+r)/2)==0
    assert s.factor(moments['omega2']-2*(1-r))==0
    prior=json.loads((ROOT/'results/susskind-memory-oriented-kernel/reverse-cross.json').read_text())
    parse=lambda z:s.sympify(z,locals={'q':q,'lam':lam,'C':C})
    V4,V2=[parse(prior[k]) for k in ['V4','V2']]
    W=4*C**2*lam+16*C**2*lam**2*q/(1+q);VF=64*lam**3*C**3*q*q
    D3=(1+q*q)*(1+q+q*q)
    A1=4*(1-q*q)*(1+q*q)*(1+q+q*q)**2
    A2=16*q*q*(1+q)*(1+q*q)*(1+q+q*q)*(2*q*q-q+2)
    VR=C**3*(lam*A1+lam**2*A2)/D3**2
    bJ,bD,bR,bF,kap=[co[k] for k in ['bJ','bD','bR','bF','kappa']]
    H=bJ*e**2+bD*om**2+kap # Dimensionless: actual energy multiplier is C*H.
    hmean=bJ*moments['e2']+bD*moments['omega2']+kap
    hsq=bJ*bJ*moments['e4']+bD*bD*moments['omega4']+kap*kap+2*bJ*bD*moments['e2omega2']+2*kap*(bJ*moments['e2']+bD*moments['omega2'])
    whsq=bJ*bJ*moments['e4omega2']+bD*bD*moments['omega6']+kap*kap*moments['omega2']+2*bJ*bD*moments['e2omega4']+2*kap*(bJ*moments['e2omega2']+bD*moments['omega4'])
    norm_parts=[V4,V2,bR*bR*VR*C*moments['e2'].subs(r,q**3),
                bF*bF*VF*C*moments['e2'].subs(r,q),W*C**2*hsq.subs(r,q*q)]
    curv_parts=[2*C*(1-q**4)*V4,2*C*(1-q*q)*V2,
                bR*bR*VR*C**2*moments['e2omega2'].subs(r,q**3),
                bF*bF*VF*C**2*moments['e2omega2'].subs(r,q),W*C**3*whsq.subs(r,q*q)]
    cross=W*C*hmean.subs(r,q*q)
    # Series each factor before products avoids expanding enormous rational forms.
    t=s.symbols('t');order=12
    seriesq=s.exp(-4*t).series(t,0,order).removeO()
    def leading(expr):
        return s.limit((lam*expr.subs(C,1/lam)).subs(q,s.exp(-4*lam)),lam,0,dir='+')
    norm_leads=[leading(v) for v in norm_parts]
    curv_leads=[leading(v) for v in curv_parts]
    cross_lead=leading(cross)
    assert norm_leads==[s.Rational(768,5),0,0,0,s.Rational(352,5)],norm_leads
    assert sum(curv_leads)==7680 and cross_lead==16
    assert s.factor(sum(curv_leads)/sum(norm_leads)-s.Rational(240,7))==0
    # Derivative identity for f_j=sech^(2j)(2 Js T), Js=1 here.
    fs=s.symbols('f1:6')
    second=lambda j:16*j*j*fs[j-1]-8*j*(2*j+1)*fs[j]
    d4f2=s.expand(64*second(2)-80*second(3))
    total=s.expand(s.Rational(768,5)*fs[3]+(d4f2-128*second(2)+4096*fs[1])/100)
    assert s.expand(total-(288*fs[3]-64*fs[2]))==0
    data=dict(status='passed',coefficients={k:str(s.factor(v)) for k,v in co.items()},
        free_algebra_raw_identity_coordinates=list(map(str,diff)),
        unit_variance_probe_moments={k:str(v) for k,v in moments.items()},
        norm_parts=list(map(str,norm_parts)),curvature_parts=list(map(str,curv_parts)),YU_equal_time=str(cross),
        ordered_small_lambda=dict(lambda_norm_parts_Js8=list(map(str,norm_leads)),
          lambda_curvature_parts_Js10=list(map(str,curv_leads)),lambda_YU_Js6=str(cross_lead),
          normalized_curvature='240/7',normalized_kernel='(9 sech^8(2 Js T)-2 sech^6(2 Js T))/7',
          projected_Q4_norm_fraction='24/35',Y_descendant_norm_fraction='11/35'),
        environment=dict(python=platform.python_version(),sympy=s.__version__))
    names=['code/memory_raw_spectrum.py','results/susskind-memory-descendants/symbolic-filters.json',
           'results/susskind-memory-oriented-kernel/reverse-cross.json']
    data['input_sha256']={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names}
    target.write_text(json.dumps(data,indent=2)+'\n')
    print('passed raw identity, spectral moments; norm leads',norm_leads,'curvature',curv_leads,flush=True)


if __name__=='__main__':main()
