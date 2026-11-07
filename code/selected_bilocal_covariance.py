"""Pure two-field selected response and one-U Gram calculus, C=1."""
from fractions import Fraction as F
from itertools import product
from functools import lru_cache
from slater_higher_wick_sources import moment
from memory_joint_regression import JointModel,Term
from neutral_tower_projection import kappa,tower,to_symbol,add
from neutral_selected_response import chain_moment

@lru_cache(None)
def term_pair(s,t,a,b,gamma,lam,mode):
    """tau/sigma(s^dag h^a t h^b), preserving both OU/chaos regions."""
    s=s.dagger();offset=s.length+a;length=offset+t.length+b
    shift=lambda items:tuple((x+offset,y+offset,*tail) for x,y,*tail in items)
    marks=s.marks+shift(t.marks);wr=s.wicks+shift(t.wicks)
    filters=s.filters+shift(t.filters)
    assert not s.xpower and not t.xpower
    ranges=[([sel] if sel>=0 else range((y-x)%2,y-x+1,2)) for x,y,p,sel in filters]
    result=F(0)
    for ds in product(*ranges):
        weight=gamma**(4*sum(p*d for (_,_,p,_),d in zip(filters,ds)))
        cr=tuple((x,y,d) for (x,y,_,_),d in zip(filters,ds))
        result+=weight*moment(length,marks,gamma,lam,wick_regions=wr,chaos_regions=cr,mode=mode)
    return result

def field_pair(left,right,a,b,gamma,lam,mode='trace'):
    return sum((u*v*term_pair(s,t,a,b,gamma,lam,mode)
                for s,u in left.items() for t,v in right.items()),F(0))

def tilted_pair(left,right,a,b,gamma,lam):
    return ((1-gamma**2)*field_pair(left,right,a,b,gamma,lam)+
            gamma**2*field_pair(left,right,a,b+2,gamma,lam))

def chain_correction(px,py,a,b,q,lam,zeta):
    return zeta*lam**2*(1+q)*sum((v*w*chain_moment(j,i+a+k,l+b,q)
                    for (i,j),v in px.items() for (k,l),w in py.items()),F(0))

def centered_prefix(n,v,model,kind='L'):
    assert kind=='L' or v==n
    raw={Term(n,((0,n),),((0,v,kind),)):F(1)}
    reg=model.scale(model.Hermite(n),-kappa(v,model.q,model.lam,kind=='D'))
    phi=to_symbol({(i,n-1-i):F(1) for i in range(v)},model.q) if kind=='L' else {}
    return model.add(raw,reg),phi

def neutral_fields(maximum,model):
    from memory_wick_descendant_chords import descendants
    _,V,_=descendants(maximum)
    _,Xi,_,reg=tower(maximum,model.q,model.lam)
    fields=[]
    for j in range(maximum+1):
        field=model.neutral(V[j])
        for d,c in reg[j].items():field=model.add(field,model.scale(model.Hermite(d),-c))
        fields.append((field,to_symbol(Xi[j],model.q)))
    return fields

@lru_cache(None)
def term_boundary(s,t,a,b,word,gamma,lam):
    """Trace with an explicit leading h/W word, independent of selected boundary."""
    s=s.dagger();start=len(word);offset=start+s.length+a
    shift=lambda items,n:tuple((x+n,y+n,*tail) for x,y,*tail in items)
    wm=tuple((i,i+1,'L') for i,c in enumerate(word) if c=='W')
    marks=wm+shift(s.marks,start)+shift(t.marks,offset)
    wr=shift(s.wicks,start)+shift(t.wicks,offset)
    filters=shift(s.filters,start)+shift(t.filters,offset)
    ranges=[([sel] if sel>=0 else range((y-x)%2,y-x+1,2)) for x,y,p,sel in filters]
    result=F(0)
    for ds in product(*ranges):
        weight=gamma**(4*sum(p*d for (_,_,p,_),d in zip(filters,ds)))
        cr=tuple((x,y,d) for (x,y,_,_),d in zip(filters,ds))
        result+=weight*moment(offset+t.length+b,marks,gamma,lam,wick_regions=wr,
                    chaos_regions=cr,shifts=(lam,)*len(wm)+(F(0),)*(len(marks)-len(wm)))
    return result

def boundary_pair(left,right,a,b,word,gamma,lam):
    return sum((u*v*term_boundary(s,t,a,b,tuple(word),gamma,lam)
                for s,u in left.items() for t,v in right.items()),F(0))
