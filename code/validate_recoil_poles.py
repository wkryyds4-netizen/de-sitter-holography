"""Independent representations, pole limits, and direct-C0 recoil calibration."""
from recoil_poles import *
from recoil_vertex import logF

def run():
    mp.mp.dps = 40
    checks=[]; identities=[]; pinch_rows=[]; c_rows=[]; shift_rows=[]
    def check(name, error, tol):
        error=abs(error); tol=mp.mpf(tol)
        checks.append(dict(name=name,error=str(error),tolerance=str(tol),passed=bool(error<=tol)))
    es=sp.symbols('e0:6'); aa,rr=parameters(es)
    deficits=[sp.expand(sum(es[k] for k in vert)/2-es[j]) for vert in VERT for j in vert]
    differences=[sp.expand(y-x) for x in aa for y in rr]
    check('twelve-triangle-pinch-bijection',int(sorted(map(str,deficits))!=sorted(map(str,differences))),0)
    check('balanced-a-r',int(sp.simplify(sum(aa)-sum(rr))),0)
    for j in range(4):
        excess=sum(1+sp.I*(aa[j]-aa[k]) for k in range(4) if k!=j)
        excess-=sp.Rational(1,2)+sp.I*aa[j]+sum(sp.Rational(1,2)+sp.I*(aa[j]-y) for y in rr)
        check('left-excess-'+str(j),int(sp.simplify(excess-1)),0)
    for l in range(3):
        excess=1-sp.I*rr[l]+sum(1+sp.I*(rr[h]-rr[l]) for h in range(3) if h!=l)
        excess-=sum(sp.Rational(1,2)+sp.I*(x-rr[l]) for x in aa)
        check('right-excess-'+str(l),int(sp.simplify(excess-1)),0)
    M=mp.mpf(3);nu=mp.mpf(1);u=mp.mpf(1)/8;v=mp.mpf(1)/4
    for index,w in enumerate([mp.mpf('.4'),mp.mpc('.4','.6'),mp.mpc('.3','1.4')]):
        a,r,edges,ts=setup(M,nu,w,u,v)
        h,cl,_=regular_left(a,r);hr,cr,_=regular_right(a,r)
        raw,rawsum,_=residue(a,r,40)
        check('left-right-'+str(index),hr/h-1,'1e-27')
        check('raw-regular-'+str(index),raw/(h*gamma_product(a,r))-1,'1e-27')
        verts=mp.exp(sum(logT(t) for t in ts))
        sine=mp.pi**-6*leg_product(edges)*mp.fprod(mp.sin(mp.pi*p) for row in ps(a,r) for p in row)
        check('vertex-sine-'+str(index),verts/sine-1,'1e-32')
        ac,rc,_,_=setup(M,nu,mp.conj(w),u,v)
        bc=mp.conj(residue(ac,rc,40)[0])
        direct=verts*raw*bc
        regular,_,_,cond=regular_graph(M,nu,w,u,v)
        check('full-graph-reflection-'+str(index),direct/regular-1,'1e-27')
        row=dict(omega=str(w),left_H=str(h),right_H=str(hr),raw_B=str(raw),
            raw_cancellation=str(rawsum/abs(raw)),left_cancellation=str(cl),right_cancellation=str(cr),
            direct_graph=str(direct),regular_graph=str(regular),paired_condition=str(cond))
        if index==0:
            bquad=vertical(a,r,40)
            check('original-vertical-versus-residues',bquad/raw-1,'1e-27')
            row['original_vertical']=str(bquad)
        identities.append(row)
        print('representation',index,'done',flush=True)

    # Dyadic masses keep p=-N exact, including multiplicities 1, 2 and 4.
    for name,uu,vv in [('simple',u,v),('double',mp.mpf(0),v),('quadruple',mp.mpf(0),mp.mpf(0))]:
        a0,_,_,_=setup(M,nu,mp.mpf(0),uu,vv)
        wj=M+nu-a0[0]
        for N in [0,1]:
            w0=2*wj+1j*(2*N+1); a,r,_,_=setup(M,nu,w0,uu,vv)
            l=2; sign=-1; p=ps(a,r)
            matching=[j for j in range(4) if p[j][l]==-N];q=len(matching)
            expected={'simple':1,'double':2,'quadruple':4}[name]
            check(f'pinch-multiplicity-{name}-{N}',q-expected,0)
            z0,h0,hs0,c0=regular_graph(M,nu,w0,uu,vv)
            leading=h0*(2j*sign*(-1)**N/mp.factorial(N))**q*mp.fprod(
                mp.gamma(p[j][ell]) for j in range(4) for ell in range(3)
                if not (ell==l and j in matching))
            if q==1:
                analytic=simple_residue(a,r,0,l,N,sign)
                check(f'closed-simple-residue-{N}',leading/analytic-1,'1e-30')
            rows=[]
            for eps in map(mp.mpf,['.001','.00001','.0000001']):
                w=w0+eps;ap,rp,_,_=setup(M,nu,w,uu,vv)
                z,h,hs,cond=regular_graph(M,nu,w,uu,vv)
                b=h*gamma_product(ap,rp)
                rows.append(dict(epsilon=str(eps),B=str(b),H=str(h),H_sharp=str(hs),graph=str(z),
                    scaled_pole_coefficient=str(b*eps**q),relative_pole_limit=str(abs(b*eps**q/leading-1)),
                    relative_graph_limit=str(abs(z/z0-1)),condition=str(cond)))
            check(f'pole-limit-{name}-{N}',mp.mpf(rows[-1]['relative_pole_limit']),'1e-4')
            check(f'graph-removable-{name}-{N}',mp.mpf(rows[-1]['relative_graph_limit']),'1e-4')
            for metric in ['relative_pole_limit','relative_graph_limit']:
                ratio=mp.mpf(rows[-1][metric])/mp.mpf(rows[-2][metric])
                check(f'limit-refinement-{metric}-{name}-{N}',ratio,'.02')
            with mp.workdps(60):
                zh,_,_,_=regular_graph(M,nu,w0,uu,vv)
            check(f'pole-value-precision-{name}-{N}',z0/zh-1,'1e-28')
            pinch_rows.append(dict(family=name,M=str(M),nu=str(nu),u=str(uu),v=str(vv),N=N,
                omega0=str(w0),multiplicity=q,H0=str(h0),H_sharp0=str(hs0),graph0=str(z0),
                graph0_high=str(zh),leading_B=str(leading),condition=str(c0),approaches=rows))
            print('pinch',name,N,'done',flush=True)

    ref=[M,M,nu]
    for index,ms in enumerate([[M+mp.mpf('.2'),M-mp.mpf('.2'),nu],
                               [mp.mpf('3.2'),mp.mpf('3.4'),mp.mpf('1.2')],
                               [M+mp.mpc('.2','.6'),M-mp.mpc('.2','.6'),nu]]):
        c=direct_c(ms);cs=mp.conj(direct_c([mp.conj(x) for x in ms]))
        pair=direct_c_pair(ms)
        check('direct-C-Schwarz-'+str(index),c*cs/pair-1,'1e-32')
        check('C-from-V-and-U-'+str(index),v2_limit(ms)/mp.exp(logF(ms))/c-1,'1e-32')
        c_rows.append(dict(masses=list(map(str,ms)),chiral_C=str(c),sharp_C=str(cs),paired_C=str(pair)))
        if index<2:
            target=v2_limit(ms)/v2_limit(ref)
            for k in [40,160,640]:
                val=v2_finite(ms,k)/v2_finite(ref,k)
                with mp.workdps(60):hi=v2_finite(ms,k)/v2_finite(ref,k)
                check(f'V2-shift-precision-{index}-{k}',val/hi-1,'1e-27')
                shift_rows.append(dict(masses=list(map(str,ms)),kappa=k,value=str(val),high=str(hi),
                    recoil_target=str(target),relative_to_limit=str(abs(val/target-1))))

    c_poles=[]
    for N in [0,1]:
        d0=nu+1j*(2*N+1);m1=M+d0/2;m2=M-d0/2;S=2*M+nu
        leading=(-2*mp.exp(mp.pi*S)*mp.sinh(mp.pi*m1)*mp.sinh(mp.pi*m2)
            /(m1*m2*nu*mp.cosh(mp.pi*S/2)*mp.cosh(mp.pi*(2*M-nu)/2)))
        rows=[]
        for eps in map(mp.mpf,['.001','.00001','.0000001']):
            d=d0+eps;val=direct_c_pair([M+d/2,M-d/2,nu])
            rows.append(dict(epsilon=str(eps),value=str(val),scaled=str(eps*val),relative=str(abs(eps*val/leading-1))))
        check('direct-C-pole-limit-'+str(N),mp.mpf(rows[-1]['relative']),'1e-4')
        c_poles.append(dict(N=N,delta0=str(d0),residue=str(leading),approaches=rows,
            scope='bare direct-C0 factor only; unknown insertion normalization not supplied'))
    inputs=['code/recoil_poles.py','code/validate_recoil_poles.py','code/recoil_internal.py','code/recoil_vertex.py',
        'sources/chandra-liska-post-2610.30247v1.html','sources/chen-stanford-tang-yang-2607.14042v1.html']
    data=dict(created_at_utc=datetime.now(timezone.utc).isoformat(),status='passed' if all(c['passed'] for c in checks) else 'failed',
        scope='recoil topology graph Gamma-pole cancellation and bare direct-C0 comparison; not physical OTOC or interval certification',
        environment=dict(python=sys.version,platform=platform.platform(),mpmath=mp.__version__,sympy=sp.__version__,dps=40,crosscheck_dps=60),
        checks=checks,identity_rows=identities,pinch_rows=pinch_rows,direct_c_rows=c_rows,shift_rows=shift_rows,direct_c_poles=c_poles,
        hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs})
    out=ROOT/'results/stage4/recoil-poles';out.mkdir(exist_ok=True,parents=True)
    dest=out/'diagnostics.json'
    assert not dest.exists(), 'Preserve previous diagnostics before a revised run.'
    dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=data['status'],checks=len(checks),failed=[c for c in checks if not c['passed']]),indent=2),flush=True)
    if data['status']!='passed':raise SystemExit(1)

if __name__=='__main__':run()
