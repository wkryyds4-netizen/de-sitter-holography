"""Independent Decimal replay of saved group-pairing algebra (no quadratures)."""
from decimal import Decimal as D, localcontext
from math import factorial
import json

PI = D("3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628")


def audit(root, number, norm):
    errors = []
    base = root / "results/stage4/ds3-group-pairing"
    main = json.loads((base / "diagnostics.json").read_text())
    off = json.loads((base / "offdiagonal.json").read_text())

    def require(ok, message):
        if not ok:
            errors.append("group-pairing: " + message)

    def close(a, b, label, tol="1e-37"):
        require(abs(a-b) <= D(tol)*max(abs(b), D("1e-60")), label)

    def ch(x):
        e = x.exp()
        return (e+1/e)/2

    def sh(x):
        e = x.exp()
        return (e-1/e)/2

    def th(x):
        e = (2*x).exp()
        return (e-1)/(e+1)

    def density(x):
        return 2*x*x/PI

    def mu(x):
        return x*th(PI*x)/PI

    def denom(a, b, p):
        return (ch(PI*(a+b))+ch(PI*p))*(ch(PI*(a-b))+ch(PI*p))/4

    def kernel(a, b, p):
        return mu(a)*mu(b)/(8*denom(a,b,p))

    def triple(a, b, p):
        return PI*sh(PI*a)*sh(PI*b)*sh(PI*p)/(8*a*b*p*denom(a,b,p))

    def state(s):
        q = (-2*PI*s).exp()
        # q <= exp(-.4*pi) in the saved samples; 240 terms exceed 90-digit needs.
        li = sum((-q)**k/D(k*k) for k in range(1,241))
        return 4*s*s/PI**2*q/(1+q)+4*s/PI**3*(1+q).ln()-2/PI**4*li

    def tail(L, s):
        beta = 2*PI
        return 2/PI**2*(-beta*L).exp()*(
            (L**3-s*s*L)/beta+(3*L*L-s*s)/beta**2+6*L/beta**3+6/beta**4)

    with localcontext() as ctx:
        ctx.prec = 90
        for data, count in [(main,81),(off,16)]:
            require(data["status"] == "passed" and len(data["checks"]) == count
                    and len({c["name"] for c in data["checks"]}) == count,
                    "count/status/duplicate checks")
            for c in data["checks"]:
                a,b = number(c["observed"]),number(c["expected"])
                scale,err,tol = map(D,(c["scale"],c["scaled_error"],c["tolerance"]))
                actual = norm((a[0]-b[0],a[1]-b[1]))/scale
                require(scale>0 and c["passed"] is True and 0<=err<tol
                        and actual<tol and abs(actual-err)<D("1e-35"),
                        "saved check replay "+c["name"])
        for key,count in [("spherical",9),("pair_packets",2),("triples",3),("gluing",15),
                          ("joint",3),("state_norms",3),("rates",3),("means",9),("time_integrals",8)]:
            require(len(main[key])==count, "raw row count "+key)
        require(len(off["rows"])==4, "offdiagonal raw row count")
        for r in main["pair_packets"]:
            a,b,R=map(D,(r["a"],r["b"],r["R"]))
            exact=PI/2*(b-a)**17*D(factorial(8)**2)/D(factorial(17))
            coeff=[24,-480,2160,-3360,1680]
            energy=sum(D(x*y)/D(i+j+1) for i,x in enumerate(coeff) for j,y in enumerate(coeff))
            bound=PI/2*(b-a)**9*energy/R**8
            close(D(r["schur"]),exact,"pair target")
            close(D(r["tail_bound"]),bound,"pair tail")
            error=abs(D(r["radial"])-exact)
            require(error<=bound, "pair integral outside tail")
            require(abs(error-D(r["error"]))<D("1e-48"),"pair saved error")
        for r in main["triples"]:
            a,b,p=map(D,r["nus"]);S=a+b+p
            exact=triple(a,b,p)
            tanh=PI/(8*a*b*p)*(sum(th(PI*(S-2*n)/2) for n in [a,b,p])-th(PI*S/2))
            close(D(r["closed"]),exact,"K3 product")
            close(D(r["tanh"]),tanh,"K3 tanh")
            close(exact,tanh,"K3 identity")
            q=D(-65).exp()
            bound=((1+q)/(1-q)).ln()/(a*b*p)
            close(D(r["tail_bound"]),bound,"K3 tail")
            require(abs(D(r["direct"])-exact)<=bound+D("1e-45"),"K3 direct integral")
        for r in main["gluing"]:
            z,w=number(r["left"]),number(r["right"])
            require(norm((z[0]-w[0],z[1]-w[1]))<D("1e-36")*max(norm(w),D("1e-60")),
                    "angular gluing")
        for r in main["joint"]:
            a,b,p=map(D,r["nus"])
            k3=triple(a,b,p);Q=kernel(a,b,p)
            f=1/ch(PI*b)
            B=PI/(4*(4*PI).sqrt()*denom(a,b,p))*(sh(PI*a)*sh(PI*b)*sh(PI*p)/(a*b*p)).sqrt()
            E=(mu(a)*mu(b)).sqrt()*B*ch(PI*b)
            close(D(r["state_amplitude"]),f,"linear state")
            close(D(r["emission"]),E,"emission amplitude")
            close(D(r["triple_norm"]),k3,"joint K3")
            close(D(r["raw_joint"]),f*f*E*E/k3,"raw quotient")
            close(D(r["closed"]),Q,"physical Q")
            M=(a+b)/2;w=a-b;A=ch(2*PI*M);B0=ch(PI*w);C=ch(PI*p)
            old=M*M*(-2*PI*M).exp()/(8*PI*(B0+C))*(1-w*w/(4*M*M))*(A-B0)/(A+C)
            n=1/((1+(-2*PI*a).exp())*(1+(-2*PI*b).exp()))
            close(D(r["old_joint"]),old,"old kernel")
            close(D(r["ratio"]),PI*Q/(8*old),"old ratio")
            close(D(r["n"]),n,"local n")
            close(PI*Q/(8*old),n,"local n identity")
        for r in main["state_norms"]:
            z=state(D(r["cutoff"]))
            close(D(r["analytic"]),z,"state analytic normalization")
            close(D(r["direct"]),z,"state quadrature normalization")
        for r in main["rates"]:
            b,s=map(D,(r["b"],r["cutoff"]))
            bound=b/(2*PI)+1/(2*PI**2)+1/(4*PI**3*s)
            close(D(r["bound"]),bound,"creation rate bound")
            require(0<D(r["rate"])<=bound,"creation rate inequality")
        for r in main["means"]:
            M,w,p=map(D,(r["M"],r["w"],r["nu"]));a=M+w/2;b=M-w/2
            A,B,C=ch(2*PI*M),ch(PI*w),ch(PI*p)
            R=(1-w*w/(4*M*M))*(A-B)/(A+C)*(A+1)/(A+B)
            reduced=kernel(a,b,p)/(density(M)/ch(PI*M)**2)
            close(D(r["reduced"]),reduced,"mean kernel")
            close(D(r["expected"]),R/(8*PI*(B+C)),"mean expected")
            close(D(r["correction"]),R,"mean correction")
            delta=(1+a*a).sqrt()-(1+b*b).sqrt();bound=abs(w)/(2*min(a,b)**2)
            require(abs(D(r["proper_difference"])-delta)<D("1e-37"),"proper energy")
            require(abs(delta-w)<=bound and 0<R<=1,"proper energy/correction bound")
            require(abs(D(r["energy_error_bound"])-bound)<D("1e-37"),"saved energy bound")
        for r in main["time_integrals"]:
            coarse,refined=number(r["coarse"]),number(r["refined"])
            require(norm((coarse[0]-refined[0],coarse[1]-refined[1]))<D("2e-10")*norm(refined),
                    "time refinement")
            require(abs(refined[1])<D("1e-14"),"time imaginary part")
            require(r["orders"]==[80,120],"time orders")
            L,s=map(lambda x:D(str(x)),(r["L"],r["cutoff"]))
            close(D(r["analytic_outer_tail"]),tail(L,s)/state(s),"time outer tail")
        for r in off["rows"]:
            a,b,p=map(D,r["nus"]);l=D(r["l"]);T=D(r["T"])
            k=triple(a,b,p);v=number(r["invariant"]);expected=tuple(k*x for x in v)
            saved=number(r["expected"]);direct=number(r["direct"])
            close(D(r["spherical_group"]),k,"offdiagonal K3")
            require(norm(tuple(x-y for x,y in zip(saved,expected)))<D("1e-37")*norm(expected),
                    "offdiagonal rank-one product")
            bound=2*(2*l+1).sqrt()/p/(1-(-2*T).exp())*(-T).exp()*(T*T+2*T+2)
            close(D(r["tail_bound"]),bound,"offdiagonal tail")
            require(norm(tuple(x-y for x,y in zip(direct,expected)))<=bound,"offdiagonal integral")
    return errors
