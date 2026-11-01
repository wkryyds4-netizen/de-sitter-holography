"""Independent Wick-sum checks and exact saved endpoint arithmetic."""
from pathlib import Path
from fractions import Fraction as F
from math import factorial,comb
import sys,json,hashlib
root=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(root/'tmp/arb-runtime'))
from flint import arb,acb,ctx

def bounds(d):return F(*map(int,d['lower'])),F(*map(int,d['upper']))

def exact_moment(m,n,sign):
    real=F(0);imag=F(0)
    if (m+n)%2:return real,imag
    for k in range(min(m,n)+1):
        if (m-k)%2 or (n-k)%2:continue
        a,b=(m-k)//2,(n-k)//2
        # Number of cross pairs times independent same-coordinate pairings.
        count=comb(m,k)*comb(n,k)*factorial(k)
        count*=factorial(2*a)//(2**a*factorial(a))
        count*=factorial(2*b)//(2**b*factorial(b))
        value=F(count*2**(a+b)*sign**k,5**((m+n)//2))
        if k%4==0:real+=value
        elif k%4==1:imag+=value
        elif k%4==2:real-=value
        else:imag-=value
    return real,imag

def audit():
    errors=[];count=0
    def ck(x,msg):
        nonlocal count
        count+=1
        if not x:errors.append(msg)
    base=root/'results/stage4/residual-certificate'
    data=json.loads((base/'chirp.json').read_text())
    ck(data['status']=='chirp steps certified; F-step overlaps uncomputed','scope')
    for p,h in data['hashes'].items():ck(hashlib.sha256((root/p).read_bytes()).hexdigest()==h,'hash '+p)
    for row in data['rows']:
        ck(hashlib.sha256((root/row['input_file']).read_bytes()).hexdigest()==row['input_sha256'],'input hash')
        for step in row['chirp_steps']:
            lo,hi=bounds(step['residual_squared']);a,b=bounds(step['residual'])
            ck(0<lo<=hi and hi-lo<F(1,10**24),'narrow positive residual squared')
            ck(a*a<=hi and b*b>=lo,'square-root interval compatibility')
        for c in row['conditional_margins']:
            ck('uncomputed' in c['hypothesis'],'conditional label')
            lo,hi=bounds(c['margin']);ck(c['positive']==(lo>0),'conditional sign')
    for a,b in zip(data['rows'][1]['chirp_steps'],data['rows'][2]['chirp_steps']):
        for key in ['overlap_real','overlap_imag','residual_squared','residual']:
            al,ah=bounds(a[key]);bl,bh=bounds(b[key])
            ck(al<=bl<=bh<=ah,'768/1024 nested interval '+key)
    initial=json.loads((base/'chirp-initial-tail-control.json').read_text())
    old=hashlib.sha256((root/'code/chirp_residual_certificate_initial.py').read_bytes()).hexdigest()
    ck(initial['hashes']['code/chirp_residual_certificate.py']==old,'historical code alias')
    for oldrow,newrow in zip(initial['rows'],data['rows']):
        for oldstep,newstep in zip(oldrow['chirp_steps'],newrow['chirp_steps']):
            ck(all(oldstep[k]==newstep[k] for k in ['overlap_real','overlap_imag','residual_squared','residual']),
               'corrected control leaves Gaussian moments unchanged')
    # Independently sum pairings instead of trusting the production recurrence.
    ctx.prec=1024
    samples=[(0,0),(2,0),(1,1),(2,2),(7,9),(8,12),(31,33),(64,64),
             (126,126),(127,127),(254,0),(0,254),(254,254),(127,128)]
    size=255
    for sign in [-1,1]:
        a=arb(2)/5;c=acb(0,arb(sign)/5)
        mom=[[acb(0) for _ in range(size)] for _ in range(size)];mom[0][0]=acb(1)
        for j in range(2,size,2):mom[0][j]=(j-1)*a*mom[0][j-2]
        for i in range(1,size):
            for j in range(size):
                mom[i][j]=((i-1)*a*mom[i-2][j] if i>1 else 0)+(j*c*mom[i-1][j-1] if j else 0)
        for m,n in samples:
            real,imag=exact_moment(m,n,sign)
            for value,ball in [(real,mom[m][n].real),(imag,mom[m][n].imag)]:
                lo=ball.lower().fmpq();hi=ball.upper().fmpq()
                ck(F(int(lo.p),int(lo.q))<=value<=F(int(hi.p),int(hi.q)),f'exact Wick sum {m},{n},{sign}')
    pilot=json.loads((base/'rectangle-pilot.json').read_text())
    ck('full-plane overlap uncomputed' in pilot['status'],'pilot scope')
    for p,h in pilot['hashes'].items():ck(hashlib.sha256((root/p).read_bytes()).hexdigest()==h,'pilot hash '+p)
    from residual_overlap_rectangle import outer_guard
    ck(outer_guard(acb(0)),'guard origin')
    ck(not outer_guard(acb(arb.pi().sqrt())),'guard actual pole row')
    ck(not outer_guard(acb(arb(-1)*arb.pi().sqrt())),'guard negative pole row')
    ck(not outer_guard(acb(arb(0,2))),'guard interval crosses row')
    return dict(status='passed' if not errors else 'failed',checks=count,errors=errors,
        scope='Wick recurrence, saved endpoint consistency and provenance; uncomputed F overlaps are not certified.')

if __name__=='__main__':
    out=root/'results/stage4/residual-certificate/chirp-audit.json'
    if out.exists():raise RuntimeError('Do not overwrite')
    data=audit();data['hashes']={'code/audit_chirp_residual.py':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    out.write_text(json.dumps(data,indent=2));print(json.dumps(data))
