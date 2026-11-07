"""Cubic Gaussian-chaos derivative singlet: marked-word DS Gram and kinetics."""
from pathlib import Path
from itertools import combinations
from fractions import Fraction
import json, platform, hashlib
import numpy as np
import sympy as s
from complex_syk_dictionary import fermions, creators, comm
from check_multimark_car import local_coeff
from singlet_coefficient import pairings
from singlet_spectral_graphs import graph

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'results/susskind-cubic-singlet'
q,z=s.symbols('q lambda', positive=True)
marks={
 'T':[((0,2),s.Integer(1)),((0,1),-s.Rational(1,2)),((1,2),-s.Rational(1,2))],
 'S':[((0,1),s.Rational(1,2)),((1,2),s.Rational(1,2))],
 'Z':[(None,s.Integer(1))]}

def moment(left,right,gap=0,end=0):
    length=6+gap+end; start=3+gap
    blocks=[set(range(3)),set(range(start,start+3))]
    total=0;count=0
    for pairing in pairings(tuple(range(length))):
        if any(a in block and b in block for a,b in pairing for block in blocks):continue
        count+=1
        for ma,ca in marks[left]:
            for mb,cb in marks[right]:
                mb=None if mb is None else tuple(x+start for x in mb)
                g=graph(length,pairing,(ma,mb))
                if ma is None and mb is None:v=1
                elif ma is None:v=z*g['beta']
                elif mb is None:v=z*g['alpha']
                else:v=z*g['delta']+z*z*g['alpha']*g['beta']
                total+=ca*cb*q**g['crossings']*v
    return s.factor(total),count

def car_checks():
    rng=np.random.default_rng(110611);records=[];identities=[]
    for n,p in [(3,1),(4,2)]:
        cs=fermions(n);ss=list(combinations(range(n),p));bs=dict(zip(ss,creators(cs,p)))
        dim=1<<n;Q=np.diag([i.bit_count() for i in range(dim)])
        def D(v):return Q@v-sum(c.conj().T@v@c for c in cs)
        def mark(a,b):return 2*p*a@b-D(a@b)
        def block(v):
            a,b,c=v; k12=mark(a,b)@c;k23=a@mark(b,c)
            k13=3*p*a@b@c-D(a@b@c)-k12-k23
            return [a@b@c,k12,k23,k13]
        # Integer Hermitian couplings, so all entries and identities are exact.
        basis=list(bs.values());raw=rng.integers(-2,3,size=(len(ss),len(ss)))+1j*rng.integers(-2,3,size=(len(ss),len(ss)))
        coupling=raw+raw.conj().T
        H=sum(coupling[i,j]*a@b.conj().T for i,a in enumerate(basis) for j,b in enumerate(basis))
        f=[[c,comm(H,c)] for c in cs]
        for row in f:row.extend([comm(H,row[1]),comm(H,comm(H,row[1]))])
        B=lambda a,b:sum(row[a].conj().T@row[b] for row in f)
        A=B(1,1);M=B(0,3);T=(M+M.conj().T)/2
        tests={'M1':B(0,1)+p*H,'M2':B(0,2)-A,
         'M3':M-(2*B(1,2)-B(2,1)), 'dagger':M-M.conj().T-3*comm(H,A),
         'rank_A':A-(2*p*H@H-D(H@H)),
         'rank_T':T-(3*p*H@H@H-D(H@H@H)-s.Rational(3,2)*(H@A+A@H))}
        for name,value in tests.items():
            assert np.all(value==0),(n,p,name)
            identities.append(dict(N=n,p=p,name=name))
        allowed=[pa for pa in pairings(tuple(range(6))) if not any((a<3)==(b<3) for a,b in pa)]
        for case in range(96):
            pa=allowed[case%len(allowed)];pi={a:b for a,b in pa}|{b:a for a,b in pa}
            sets=[ss[int(rng.integers(len(ss)))]]*6 if case<24 else [ss[int(rng.integers(len(ss)))] for _ in range(6)]
            ops=[bs[sets[a]]@bs[sets[pi[a]]].conj().T for a in range(6)]
            left,right=block(ops[:3]),block(ops[3:])
            for i,ma in enumerate([None,(0,1),(1,2),(0,2)]):
                for j,mb in enumerate([None,(3,4),(4,5),(3,5)]):
                    ms=[m for m in [ma,mb] if m is not None]
                    expected=local_coeff(n,p,sets,pa,ms)[-1]
                    actual=np.trace(left[i]@right[j]);assert actual.imag==0 and actual.real==round(actual.real)
                    actual=Fraction(int(actual.real),dim)
                    assert actual==expected,(n,p,case,ma,mb,actual,expected)
            records.append(dict(N=n,p=p,pairing=pa,sets=sets))
    return dict(seed=110611,operator_identities=identities,words=len(records),mark_subset_checks=16*len(records),records=records)

def main():
    OUT.mkdir(exist_ok=True)
    target=OUT/'study.json'
    if target.exists():raise FileExistsError(target)
    car=car_checks();print('CAR checks passed',flush=True)
    names=['T','S','Z'];G=s.zeros(3);K=s.zeros(3)
    for i,a in enumerate(names):
        for j,b in enumerate(names):
            G[i,j]=moment(a,b)[0]
            K[i,j]=s.factor(moment(a,b,0,2)[0]+moment(b,a,0,2)[0]-2*moment(a,b,1,1)[0])
    assert G==G.T and K==K.T
    coeff=s.simplify(G[1:,1:].inv()*G[1:,0]);v=s.Matrix([1,-coeff[0],-coeff[1]])
    variance=s.factor((v.T*G*v)[0]);kinetic=s.factor((v.T*K*v)[0])
    mu2=s.factor(kinetic/variance/z) # Js=1, sigma^2=1/lambda
    print('G=',G,'\nK=',K,'\ncoeff=',coeff,'\nvariance=',variance,'\nkinetic=',kinetic,'\nmu2=',mu2,flush=True)
    data=dict(status='passed',conventions='full trace beta=mu=0; DS at fixed lambda; q=exp(-4lambda); G/sigma^6,K/sigma^8; sigma^2=Js^2/lambda',
      car=car,gram=[[str(G[i,j]) for j in range(3)] for i in range(3)],
      kinetic=[[str(K[i,j]) for j in range(3)] for i in range(3)],
      projection=list(map(str,coeff)),residual_variance=str(variance),residual_kinetic=str(kinetic),mu2_over_Js2=str(mu2),
      small_lambda_variance=str(s.series(variance.subs(q,s.exp(-4*z)),z,0,4)),
      small_lambda_kinetic=str(s.series(kinetic.subs(q,s.exp(-4*z)),z,0,5)),
      small_lambda_mu2=str(s.limit(mu2.subs(q,s.exp(-4*z)),z,0)),
      pairing_counts=dict(gram=moment('Z','Z')[1],adjacent_H=moment('Z','Z',0,2)[1],split_H=moment('Z','Z',1,1)[1]),
      environment=dict(python=platform.python_version(),numpy=np.__version__,sympy=s.__version__))
    inputs=['code/cubic_singlet.py','code/check_multimark_car.py','code/singlet_coefficient.py','code/singlet_spectral_graphs.py','code/complex_syk_dictionary.py','code/singlet_multimark.py']
    data['input_sha256']={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs}
    target.write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps({k:data[k] for k in ['small_lambda_variance','small_lambda_kinetic','small_lambda_mu2','pairing_counts']}),flush=True)

if __name__=='__main__':main()
