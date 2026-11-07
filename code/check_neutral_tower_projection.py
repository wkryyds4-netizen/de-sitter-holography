"""Independent finite checks of projection, chaos filters and actual descendants."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
import hashlib,json,platform,time
from slater_higher_wick_sources import moment
from slater_joint_site_polynomial import pairings,crossed
from memory_wick_descendant_chords import descendants
from memory_joint_regression import JointModel
from neutral_tower_projection import *
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/susskind-selected-bilocal-mean/projection-checks.json'
if OUT.exists():raise FileExistsError(OUT)

def colored(w,wick,q,lam):
    ans=F(0)
    for pi in pairings(tuple(range(len(w)))):
        if any(w[i]!=w[j] or any(a<=i<j<b for a,b in wick) for i,j in pi):continue
        cr=sum(crossed(e,f) for i,e in enumerate(pi) for f in pi[i+1:])
        ans+=q**cr*lam**sum(w[i]=='W' for i,j in pi)
    return ans

def external(word,offset):
    return tuple((offset+i,offset+i+1,'L') for i,c in enumerate(word) if c=='W')

def source_word(n,v,test,gamma,lam,kind='L',prefix=0,suffix=0,chaos=-1,ou=0):
    length=prefix+n+suffix
    marks=((prefix,prefix+v,kind),)+external(test,length)
    # A filter over the whole field is a finite exact sum of chaos selections.
    ds=[chaos] if chaos>=0 else range(length%2,length+1,2)
    result=F(0)
    for d in ds:
        result+=gamma**(4*ou*d)*moment(length+len(test),marks,gamma,lam,
            wick_regions=((prefix,prefix+n),),chaos_regions=((0,length,d),),
            shifts=(F(0),)+(lam,)*(len(marks)-1))
    return result

def projected_word(poly,test,q,lam):
    return sum((c*colored(('h',)*a+('W',)+('h',)*b+test,(),q,lam) for (a,b),c in poly.items()),F(0))

def field_moment(field,test,model,gamma,mode='trace'):
    ans=F(0)
    for term,c in model.neutral(field).items():
        assert not term.xpower
        marks=term.marks+external(test,term.length)
        shifts=(F(0),)*len(term.marks)+(model.lam,)*(len(marks)-len(term.marks))
        ds=[None] if not term.filters else range(term.length%2,term.length+1,2)
        for d in ds:
            weight=F(1);cr=[]
            if d is not None:
                for a,b,power,selected in term.filters:
                    assert (a,b)==(0,term.length)
                    if selected>=0 and d!=selected:weight=0
                    weight*=model.q**(power*d);cr.append((a,b,d))
            if weight:
                ans+=c*weight*moment(term.length+len(test),marks,gamma,model.lam,
                    wick_regions=term.wicks,chaos_regions=tuple(cr),shifts=shifts,mode=mode)
    return ans

start=time.monotonic();records=[];counts={}
def check(kind,actual,predicted,**data):
    row=dict(kind=kind,actual=str(actual),predicted=str(predicted),passed=actual==predicted,**data)
    records.append(row);counts[kind]=counts.get(kind,0)+1
    if actual!=predicted:
        OUT.parent.mkdir(exist_ok=True,parents=True)
        OUT.with_name('projection-failure.json').write_text(json.dumps(records,indent=2)+'\n')
        raise AssertionError(row)

for gamma,lam in [(F(1,2),F(1,3)),(F(2,3),F(1,5))]:
 q=gamma**4;model=JointModel(q,lam);E,Xi,g,vreg=tower(5,q,lam)
 # All two-colour words at the minimal nonzero external length, n through four.
 for n in range(1,5):
  for v in range(1,n+1):
   for test in product('hW',repeat=n):
    act=source_word(n,v,test,gamma,lam)
    pred=kappa(v,q,lam)*colored(('h',)*n+test,((0,n),),q,lam)
    pred+=sum((colored(('h',)*i+('W',)+('h',)*(n-i-1)+test,((0,n),),q,lam) for i in range(v)),F(0))
    check('prefix',act,pred,n=n,v=v,test=''.join(test),q=str(q))
   # Charged pure source has only the known energy projection.
   for test in [('h',)*n,('W',)*n]:
    check('charged',source_word(n,n,test,gamma,lam,kind='D'),
          kappa(n,q,lam,True)*colored(('h',)*n+test,((0,n),),q,lam),n=n,test=''.join(test),q=str(q))
 print(json.dumps({'q':str(q),'phase':'prefix','checks':len(records)}),flush=True)
 # Prefix multiplication creates lower chaos. Filter it before projection.
 for n,v,l,r in [(1,1,1,0),(1,1,1,1),(2,1,1,0),(2,2,1,0),(2,1,1,1)]:
  module={(i,n-i-1):F(1) for i in range(v)}
  for _ in range(l):module=left(module,q)
  for _ in range(r):module={(b,a):c for (a,b),c in left({(b,a):c for (a,b),c in module.items()},q).items()}
  length=n+l+r
  for d in range(length%2,length+1,2):
   projected=to_symbol({ab:c for ab,c in module.items() if sum(ab)+1==d},q)
   for test in [('W',)+('h',)*k for k in range(4) if (length+k+1)%2==0]:
    # Odd W test kills the energy-only component.
    check('chaos',source_word(n,v,test,gamma,lam,prefix=l,suffix=r,chaos=d),
          projected_word(projected,test,q,lam),n=n,v=v,left=l,right=r,d=d,test=''.join(test),q=str(q))
 # Original ordered-field recursion + local source evaluator versus one-U recursion.
 eta,V,_=descendants(5)
 for j in range(1,6):
  for k in range(3):
   test=('W',)+('h',)*k
   if (j+k+1)%2 or j+k+1>8:continue
   for name,field,module in [('eta',eta[j],E[j]),('V',V[j],Xi[j])]:
    check('descendant',field_moment(field,test,model,gamma),
          projected_word(to_symbol(module,q),test,q,lam),field=name,j=j,test=''.join(test),q=str(q))
 # Selected expectations keep the full density boundary, not the trace shortcut.
 for j in range(1,5):
  for a,b in [(0,0),(1,0),(0,1)]:
   if (j+a+b)%2:continue
   # Insert external h through the original field API.
   from memory_wick_descendant_chords import hleft,hright
   field=V[j]
   for _ in range(a):field=hleft(field) if not any(t.ou_power for t in field.terms) else None
   # OU terms cannot be extended by that API; test only undressed V for j>=3.
   if field is None or (b and any(t.ou_power for t in field.terms)):continue
   for _ in range(b):field=hright(field)
   energy=energy_polynomial(vreg[j],q)
   reg=sum((c*moment(k+a+b,(),gamma,lam,mode='selected') for k,c in energy.items()),F(0))
   poly={(x+a,y+b):c for (x,y),c in to_symbol(Xi[j],q).items()}
   check('selected',field_moment(field,(),model,gamma,mode='selected')-reg,
          selected_symbol(poly,q,lam,gamma),j=j,a=a,b=b,q=str(q))
 print(json.dumps({'q':str(q),'phase':'complete','counts':counts}),flush=True)

inputs=['code/neutral_tower_projection.py','code/check_neutral_tower_projection.py',
        'code/slater_higher_wick_sources.py','code/slater_composite_sources.py',
        'code/memory_wick_descendant_chords.py','code/memory_joint_regression.py',
        'code/neutral_selected_response.py']
result=dict(status='passed',scope='Exact finite polynomial/source checks; all-degree and time convergence rely on the accompanying proof.',counts=counts,rows=records,python=platform.python_version(),parameters='Two independent rational (q,lambda) pairs, not the physical exponential curve.',seconds=time.monotonic()-start,input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs})
OUT.parent.mkdir(parents=True,exist_ok=True)
with OUT.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps({'status':'passed','counts':counts,'seconds':result['seconds'],'sha256':hashlib.sha256(OUT.read_bytes()).hexdigest()}),flush=True)
