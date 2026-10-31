"""Independent Decimal replay of two-matter saved algebra and error claims."""
from decimal import Decimal as D,localcontext
from itertools import product
import json
PI=D("3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628")

def audit(root,number,norm):
    base=root/"results/stage4/ds3-two-matter"
    main=json.loads((base/"diagnostics.json").read_text())
    angular=json.loads((base/"angular.json").read_text())
    errors=[]
    def req(ok,msg):
        if not ok:errors.append("two-matter: "+msg)
    def close(a,b,msg,tol="1e-36"):
        req(abs(a-b)<=D(tol)*max(abs(b),D("1e-60")),msg)
    def ch(x):
        e=x.exp();return (e+1/e)/2
    def sh(x):
        e=x.exp();return (e-1/e)/2
    def sincos(x,cos=False):
        term=D(1) if cos else x;total=term
        for n in range(1,200):
            term*=-x*x/D((2*n-1 if cos else 2*n)*(2*n if cos else 2*n+1))
            total+=term
            if abs(term)<D("1e-85"):break
        return total
    def denom(a,b,p):return (ch(PI*(a+b))+ch(PI*p))*(ch(PI*(a-b))+ch(PI*p))/4
    def K(a,b,p):return PI*sh(PI*a)*sh(PI*b)*sh(PI*p)/(8*a*b*p*denom(a,b,p))
    def rho(b):return 2*b*b/PI
    def sigma(b):return rho(b)/ch(PI*b)**2
    def Q(a,b,p):return a*b*sh(PI*a)*sh(PI*b)/(8*PI**2*ch(PI*a)*ch(PI*b)*denom(a,b,p))
    def g(x):return x*(PI/((1+x*x).sqrt()+x)).exp()/(1+(-2*PI*x).exp())
    def bal(a,b,p):return (ch(PI*(a+b))-ch(PI*(a-b)))/(32*PI*denom(a,b,p))
    def lr(t,p):return (p if not t else sincos(p*t)/sh(t))/(4*PI*sh(PI*p))
    def deriv(t,p):
        return D(0) if not t else (p*sincos(p*t,True)/sh(t)-sincos(p*t)*ch(t)/sh(t)**2)/(4*PI*sh(PI*p))
    def tail(R,s):
        B=2*(PI*((1+s*s).sqrt()-s)).exp()
        return (-PI*R).exp()/(4*PI**2)*((1+B)*(R/PI+1/PI**2)+B/s*(R*R/PI+2*R/PI**2+2/PI**3))
    with localcontext() as ctx:
        ctx.prec=90
        for data,count in [(main,56),(angular,10)]:
            req(data["status"]=="passed" and len(data["checks"])==count and
                len({c["name"] for c in data["checks"]})==count,"check count/status")
            for c in data["checks"]:
                z,w=number(c["observed"]),number(c["expected"])
                scale,err,tol=map(D,(c["scale"],c["scaled_error"],c["tolerance"]))
                actual=norm((z[0]-w[0],z[1]-w[1]))/scale
                req(scale>0 and c["passed"] is True and 0<=err<tol and actual<tol
                    and abs(actual-err)<D("1e-35"),"saved check "+c["name"])
        for key,count in [("weights",3),("products",6),("four_spherical",3),("chains",2),("shifted",3),("limits",12)]:
            req(len(main[key])==count,"raw count "+key)
        for r in main["weights"]:
            a,q,b=map(D,r["nus"]);w=rho(b)*K(a,q,b)
            close(D(r["group"]),w,"pair group measure")
            close(D(r["ds"]),w,"DS-normalized measure")
        for r in main["products"]:
            a,q=map(D,r["nus"]);rad=D(r["r"]);L=D(r["L"])
            target=D(1) if not rad else sincos(a*rad)*sincos(q*rad)/(a*q*sh(rad)**2)
            close(D(r["target"]),target,"spherical product target")
            bound=2*sh(PI*a)/a*sh(PI*q)/q*(-PI*L).exp()*(L/PI+1/PI**2)
            close(D(r["tail"]),bound,"pair spectral tail")
            req(abs(D(r["direct"])-target)<=bound+D("1e-46"),"pair spectral integral")
        for r in main["four_spherical"]:
            nus=list(map(D,r["nus"]));total=D(0);prod=D(1)
            for n in nus:prod*=n
            for signs in product([-1,1],repeat=4):
                x=PI*sum(s*n for s,n in zip(signs,nus))/2
                f=-x*x/3 if abs(x)<D("1e-20") else 1-x*ch(x)/sh(x)
                sign=signs[0]*signs[1]*signs[2]*signs[3];total+=sign*f
            exact=total/(16*prod);close(D(r["closed"]),exact,"four spherical closed")
            T,L=map(D,(r["T"],r["L"]))
            rt=2/(prod*((2*T).exp()-1));st=(-2*PI*L).exp()
            for n in nus:st*=sh(PI*n)/n
            close(D(r["radial_tail"]),rt,"four radial tail")
            close(D(r["spectral_tail"]),st,"four spectral tail")
            req(abs(D(r["radial"])-exact)<=rt+D("1e-45"),"four radial integral")
            req(abs(D(r["spectral"])-exact)<=st+D("1e-45"),"four spectral integral")
        for r in main["chains"]:
            a,b,c,p,q=map(D,r["nus"])
            W=rho(b)*K(a,b,q)*K(b,c,p);expected=Q(a,b,q)*Q(b,c,p)/sigma(b)
            close(D(r["rigging_weight"]),W,"four rigging coefficient")
            close(norm(number(r["amplitude"]))**2/W,expected,"four amplitude quotient")
            close(D(r["expected"]),expected,"four chord kernel")
            close(D(r["raw"]),expected,"four raw norm")
        for r in main["shifted"]:
            a,b,p=map(D,r["nus"]);energy=(1+a*a).sqrt()-(1+b*b).sqrt()
            ratio=g(a)/g(b);balanced=bal(a,b,p);value=ratio*balanced
            close(D(r["energy"]),energy,"proper energy")
            close(D(r["g_ratio"]),ratio,"shift telescope factor")
            close(D(r["balanced"]),balanced,"balanced kernel")
            close(D(r["raw"]),Q(a,b,p)/sigma(b)*(PI*energy).exp(),"raw tilted kernel")
            close(D(r["tilted"]),value,"tilted identity")
            req(0<balanced<=D(r["limiting_spectrum"]),"balanced inequality")
        for r in main["limits"]:
            M,s,p,q,tp,tq,cut=map(lambda k:D(str(r[k])),["M","cutoff","p","q","tp","tq","D"])
            lead=lr(tp,p)*lr(tq,q);first=(deriv(tp,p)*lr(tq,q)+lr(tp,p)*deriv(tq,q))/M
            close(D(r["leading"]),lead,"proper four-point limit")
            z=number(r["first"])
            req(abs(z[0])<D("1e-45") and abs(z[1]-first)<D("1e-40"),"first phase correction")
            val=number(r["refined"]);coarse=number(r["coarse"])
            err=norm((val[0]-lead,val[1]));corrected=norm((val[0]-lead,val[1]-first))
            close(D(r["error"]),err,"finite-mass error","1e-35")
            close(D(r["relative_error"]),err/abs(lead),"relative error","1e-35")
            close(D(r["corrected_error"]),corrected,"first-corrected residual","1e-35")
            d=2/M+2*PI/M**2+2*PI*(-PI*M).exp();C=(M*d/2).exp()
            eps=2*(2+ch(PI*p)+ch(PI*q))*(-PI*M).exp()
            bound=C*(d+2*max(abs(tp),abs(tq))/M**2)/(2*PI**5)+eps*lr(0,p)*lr(0,q)+tail(M/2,s)
            close(D(r["uniform_absolute_bound"]),bound,"limit envelope")
            close(D(r["outer_tail"]),tail(cut,s),"difference-box tail")
            req(err<bound and r["orders"]==[64,96],"limit bound and orders")
            req(norm((val[0]-coarse[0],val[1]-coarse[1]))<D("2e-10")*norm(val),"time refinement")
        req(len(angular["rows"])==2 and len(angular["vertices"])==6,"angular counts")
        for r in angular["vertices"]:
            a,q,c,p=map(D,r["nus"]);b=D(r["b"])
            v=D(3).sqrt()*(1+a*a+q*q-b*b)/(2*((1+a*a)*(1+q*q)).sqrt())
            close(D(r["polynomial"]),v,"angular polynomial")
            z=number(r["vertex"]);req(norm((z[0]-v,z[1]))<D("1e-35"),"actual l1 vertex")
        for r in angular["rows"]:
            a,q,c,p=map(D,r["nus"]);T,L=map(D,(r["T"],r["L"]));beta=2*PI
            rt=2*D(3).sqrt()/(c*p*(1-(-2*T).exp())**2)*(-2*T).exp()*(T*T+T+D(".5"))
            st=2*PI*D(3).sqrt()/(2*((1+a*a)*(1+q*q)).sqrt())*(-beta*L).exp()
            for n in [a,q,c,p]:st*=sh(PI*n)/n
            st*=((L*L+1+a*a+q*q)/beta+2*L/beta**2+2/beta**3)
            close(D(r["radial_tail"]),rt,"angular radial tail")
            close(D(r["spectral_tail"]),st,"angular spectral tail")
            z,w=number(r["direct"]),number(r["spectral"])
            req(norm((z[0]-w[0],z[1]-w[1]))<=rt+st+D("1e-43"),"angular integral")
    return errors
