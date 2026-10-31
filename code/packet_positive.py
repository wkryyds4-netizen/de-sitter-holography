"""Independent inexpensive positive moments and direct tails of the same packet."""
import numpy as np
from packet_quadrature import inner_rule,packet_density
from quantized_observer import right_kernel,mass


def positive_at_nodes(ds,order=64,cap=12.,nu=.6,lower=.2,upper=3.,sigma=1.):
    ds=np.asarray(ds,float)
    ac,wac=inner_rule(order,lower,upper)
    # Split the tail instead of subtracting nearly equal norms.
    bounds=sorted(set([upper,min(6.,cap),cap]))
    pieces=[inner_rule(order,l,h) for l,h in zip(bounds[:-1],bounds[1:]) if h>l]
    at=np.concatenate([x for x,w in pieces]);wat=np.concatenate([w for x,w in pieces])
    def kernel(y,x):
        return right_kernel(nu,y,x)*np.exp(-(sigma*(mass(y)-mass(x)))**2)
    all_a=np.concatenate([ac,at])
    hc=wac @ kernel(ac[:,None],all_a[None,:])
    ht=wat @ kernel(at[:,None],all_a[None,:])
    first_c=kernel(ac[:,None],ds[None,:])
    first_t=kernel(at[:,None],ds[None,:])
    c2=wac @ first_c
    nc=(wac*hc[:len(ac)]) @ first_c
    inner_tail=(wac*ht[:len(ac)]) @ first_c
    outer_tail=(wat*(hc[len(ac):]+ht[len(ac):])) @ first_t
    eps=inner_tail+outer_tail
    return dict(C2_cut=c2,D_cut=c2*c2,N_cut=nc,epsilon_to_cap=eps,
        N_to_cap=nc+eps,inner_tail=inner_tail,outer_tail=outer_tail)


def average_positive(ds,weights,order=64,cap=12.):
    weights=np.asarray(weights,float)
    fields=positive_at_nodes(ds,order=order,cap=cap)
    out={key:float(weights @ value) for key,value in fields.items()}
    out["external_quadrature_mass"]=float(sum(weights))
    out["conditional_two_point_variance"]=out["D_cut"]-out["C2_cut"]**2
    out["cross_tail_bound_estimate"]=float(np.sqrt(out["epsilon_to_cap"])*(np.sqrt(out["N_to_cap"])+np.sqrt(out["N_cut"])))
    out["scope"]="Direct positive finite-cap estimates, no certified infinite remainder"
    return out


def positive_packet(outer_order=32,inner_order=64,cap=12.):
    ds,weights=inner_rule(outer_order,.7,.9)
    return average_positive(ds,weights*packet_density(ds),order=inner_order,cap=cap)