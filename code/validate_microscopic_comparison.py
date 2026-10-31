"""Check conditional coefficient identities and a pure-recoil full-function identity.

This does NOT establish equality with a microscopic four-point function.
"""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,sys
import mpmath as mp
from conformal_shock import j_boundary

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results'/'stage4'/'microscopic-comparison.json'


def main():
    mp.mp.dps=60
    checks=[];coefficients=[];full_functions=[]
    def number(x):return mp.nstr(x,48)
    def record(name,error,tolerance='1e-45'):
        checks.append({'name':name,'error':number(error),'tolerance':tolerance,
                       'passed':bool(error<mp.mpf(tolerance))})
    T=mp.mpf('.37');L=1/(2*mp.pi*T)
    def antipodal_propagator(d):
        center=mp.mpf(d-1)/2
        return mp.gamma(center+mp.mpf('.5'))*mp.gamma(center-mp.mpf('.5'))/(L**(d-2)*(4*mp.pi)**(mp.mpf(d)/2)*mp.gamma(mp.mpf(d)/2))
    for n in range(1,9):
        ratio=antipodal_propagator(3+2*n)/antipodal_propagator(3)
        # Direct d=3 substitution into Narovlansky (6.10), with the time phase
        # canceled by (-1)^n, C0 divided out, and s0^(-2n) factored out.
        microscopic=mp.gamma(n+1)**4*mp.gamma(mp.mpf('1.5')+n)**2/mp.gamma(2+2*n)**2/mp.gamma(mp.mpf('1.5'))**2
        microscopic*=ratio**2/(mp.pi**(2*n)*T**(4*n)*mp.factorial(n)**2*2**(2*n))
        # Original conformal J kernel plus independently integrated angular moment,
        # after the purely formal rescaling m=2*s0.
        moment=mp.quad(lambda phi:mp.cos(phi)**(2*n),[0,mp.pi/2,mp.pi])/mp.pi
        gravity=mp.factorial(2*n)*moment/16**(2*n)
        closed=mp.rf(mp.mpf('.5'),n)**2/64**n
        record('Eq6.10 Gamma expression n='+str(n),abs(microscopic-closed))
        record('pure-recoil angular moment n='+str(n),abs(gravity-closed))
        coefficients.append({'n':n,'microscopic_leading_log_coefficient':number(microscopic),
                             'gravity_after_formal_m_2s0':number(gravity),
                             'closed_coefficient':number(closed)})
    for av in ['0.2','0.5','1','2','10','100']:
        a=mp.mpf(av)
        plus=mp.quad(lambda phi:j_boundary(mp.cos(phi)/a,1),[0,mp.pi/2,mp.pi])/mp.pi
        minus=mp.quad(lambda phi:j_boundary(mp.cos(phi)/a,-1),[0,mp.pi/2,mp.pi])/mp.pi
        real_integral=a*mp.quad(lambda theta:mp.exp(-a*mp.sin(theta)),[0,mp.pi/2])
        residue=a*mp.besselk(0,a)
        alternative=real_integral-1j*residue
        record('E1 angle vs real integral and K0 a='+av,abs(plus-alternative))
        record('opposite pure-recoil contour a='+av,abs(minus-mp.conj(alternative)))
        full_functions.append({'a':av,'E1_upper_real':number(mp.re(plus)),'E1_upper_imag':number(mp.im(plus)),
                               'independent_real':number(real_integral),'minus_a_K0':number(-residue),
                               'absolute_error':number(abs(plus-alternative))})
    failed=sum(not r['passed'] for r in checks)
    names=['validate_microscopic_comparison.py','conformal_shock.py']
    result={'generated_utc':datetime.now(timezone.utc).isoformat(),
            'status':'failed' if failed else 'passed','checks':len(checks),'failed':failed,
            'check_details':checks,'precision_dps':60,'T_for_dimensional_check':'0.37',
            'input_sources':[
              {'url':'https://arxiv.org/html/2506.02109v2','equations':['3.6','6.10','C.39','C.40']},
              {'url':'https://arxiv.org/html/2607.14042v1','equations':['G.2','2.12','2.37','2.38','B.17']}],
            'coefficients':coefficients,'pure_recoil_full_functions':full_functions,
            'scope':'Leading-log antipodal coefficient identity after formal m=2*s0, and a pure-recoil gravity identity; NOT a microscopic full-function or observer-mass dictionary.',
            'environment':{'python':sys.version,'mpmath':mp.__version__},
            'code_sha256':{n:hashlib.sha256((ROOT/'code'/n).read_bytes()).hexdigest() for n in names}}
    OUT.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(f'{len(checks)-failed} passed; {failed} failed')
    print('maximum identity error',max(mp.mpf(c['error']) for c in checks))
    if failed:
        for c in checks:
            if not c['passed']:print(c)
        raise SystemExit(1)


if __name__=='__main__':main()
