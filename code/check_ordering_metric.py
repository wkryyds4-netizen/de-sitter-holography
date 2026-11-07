"""Independent integral and finite-q Gram controls for the limiting odd metric."""
from pathlib import Path
import json,hashlib,platform
import mpmath as mp
from two_matter_chords import TwoMatter
ROOT=Path(__file__).resolve().parent.parent;OUT=ROOT/'results/susskind-ordering-metric'

def integral_case(d,e,times,prec):
    mp.mp.dps=prec;d=mp.mpf(d);e=mp.mpf(e);u,v,up,vp=map(mp.mpf,times);D=d+e
    a,b,ap,bp=[mp.exp(2*z) for z in [u,v,up,vp]]
    c=lambda x:(a+ap)*x+(b+bp)*(1-x)
    w=lambda x:x**(2*d-1)*(1-x)**(2*e-1)
    factor=2**(2*D)*mp.exp(2*d*(u+up)+2*e*(v+vp))/(mp.gamma(2*d)*mp.gamma(2*e))
    base=mp.gamma(2*D)*mp.quad(lambda x:w(x)/c(x)**(2*D),[0,mp.mpf('.5'),1])
    deriv=mp.gamma(2*D+2)*(a-b)*(ap-bp)*mp.quad(lambda x:w(x)*x*(1-x)/c(x)**(2*D+2),[0,mp.mpf('.5'),1])
    val=factor*(deriv+2*d*e*base)
    expected=2*d*e*mp.sech(u-up)**(2*d+1)*mp.sech(v-vp)**(2*e+1)*(mp.cosh(u+up-v-vp)+mp.sinh(u-up)*mp.sinh(v-vp))
    assert abs(val-expected)<mp.mpf(10)**(-prec+12)
    return dict(d=str(d),e=str(e),times=times,dps=prec,quadratic_form_integral=str(val),kernel=str(expected),absolute_error=str(abs(val-expected)))

def finite_gram(h,prec):
    mp.mp.dps=prec;A=TwoMatter(h,2,3);B=TwoMatter(h,3,2);rows=[]
    for k in range(5):
        left=A.primary(k);right=B.primary(k);overlap=(left.T*A.matrix(k,True)*right)[0]
        coeff=2*(1-(-1)**k*overlap)/mp.mpf(h);lead=mp.mpf(k*(k+9)+12)
        rows.append(dict(k=k,overlap=str(overlap),odd_norm_over_h=str(coeff),leading=str(lead),
            relative_leading_error=str(coeff/lead-1)))
    return dict(h=h,dps=prec,records=rows)
def main():
    target=OUT/'checks.json'
    if target.exists():raise FileExistsError(target)
    integrals=[];refinements=[]
    cases=[('2','3',['0','0','0','0']),('2','3',['.4','-.2','.1','.3']),
       ('2','3',['1','0','1','0']),('2','3',['.8','-.5','-.3','.2']),
       ('1.3','2.1',['.3','-.1','-.2','.4'])]
    for d,e,t in cases:
        a=integral_case(d,e,t,55);b=integral_case(d,e,t,85)
        diff=abs(mp.mpf(a['quadratic_form_integral'])-mp.mpf(b['quadratic_form_integral']))
        assert diff<mp.mpf('1e-48');integrals.append(b);refinements.append(dict(kind='integral',times=t,difference=str(diff)))
    gram=[]
    for h in ['.04','.02','.01']:
        a=finite_gram(h,85);b=finite_gram(h,115)
        diff=max(abs(mp.mpf(x['odd_norm_over_h'])/mp.mpf(y['odd_norm_over_h'])-1) for x,y in zip(a['records'],b['records']))
        assert diff<mp.mpf('1e-65');gram.append(b);refinements.append(dict(kind='finite_gram',h=h,max_relative_change=str(diff)))
        print('gram',h,'max leading error',b['records'][-1]['relative_leading_error'],flush=True)
    mp.mp.dps=70;phase=[]
    def f(x):return (1+1j*(10*x-4)/mp.sqrt(mp.mpf(24)/11))/mp.sqrt(2)
    fp=1j*10/mp.sqrt(mp.mpf(48)/11)
    for ss in ['-.1','-.05','.05','.1']:
        st=mp.mpf(ss)
        def integrand(x):
            mprime=-(1-2*x)/(x*x*(1-x)**2)
            return 504*x**3*(1-x)**5*(x*(1-x)*abs(fp-1j*st*mprime*f(x))**2+12*abs(f(x))**2)
        got=mp.quad(integrand,[0,mp.mpf('.5'),1]);exact=17+mp.sqrt(66)*st+181*st*st
        assert abs(got-exact)<mp.mpf('1e-60');phase.append(dict(parameter=ss,norm_squared=str(got),exact=str(exact),difference=str(abs(got-exact))))
    inputs=['code/check_ordering_metric.py','code/two_matter_chords.py','code/ordering_metric.py']
    result=dict(status='passed',integrals=integrals,finite_gram=gram,phase_norm=phase,precision_refinements=refinements,
        error_scope='Finite-q Gram k<=4 and numerical integration are empirical precision checks. All-k leading operator identity is derived from the sourced covariance, not inferred from finite q overlap tests.',
        environment=dict(python=platform.python_version(),mpmath=mp.__version__),
        input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs})
    target.write_text(json.dumps(result,indent=2)+'\n');print('5 integral cases, 15 Gram cases, 4 phase norms saved',flush=True)
if __name__=='__main__':main()
