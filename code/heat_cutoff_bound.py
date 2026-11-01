"""Explicit h=1 heat-tail bound by punctured H1 approximation.

This bounds only D-D_epsilon for the specified Gaussian matrix element.
It does not bound quadrature of D_epsilon or test the pentagon itself.
"""
from pathlib import Path
import json,hashlib,platform
import mpmath as m

def evaluate(dps,epsilon,radius):
    iv=m.iv;iv.dps=dps;pi=iv.pi;rho=iv.mpf(radius);eps=iv.mpf(epsilon)
    # Seven distinct artificial product singularities near [-9,9]^2.
    poles=[(1,-3),(1,-1),(1,1),(1,3),(-3,1),(-1,1),(3,1)]
    # Gradient sums contain (pi,pi) twice, one from each product orientation.
    gradient_poles=poles+[(1,1)]
    norms=[pi*iv.sqrt(i*i+j*j) for i,j in poles]
    K=lambda d,R:iv.exp(-(d-R)**2/pi)/(pi*pi)
    mass_holes=sum(4*pi*rho*rho*K(d,2*rho) for d in norms)
    # Gaussian union tail outside the square [-8,8]^2.
    tail=iv.mpf(2)/8*iv.exp(-64/pi)
    cutoff_difference=iv.sqrt(tail+mass_holes)
    # Uniform sum of all compact-product factors j>=1 on [-9,9]^2.
    z=iv.exp(9-3*pi)
    remainder=z/((1-z)*(1-iv.exp(-2*pi)))
    constant=9*iv.sqrt(2)/pi+4+4*remainder
    inverse_distance=[]
    for i,j in gradient_poles:
        d=pi*iv.sqrt(i*i+j*j);R=iv.mpf('.5');c=iv.mpf('.6')
        near=2*pi*K(d,R)*iv.log(R/rho)
        # Split the far region at |s|=.6|p|. Gaussian radial tail is exact.
        far=1/((1-c)*d)**2+iv.exp(-(c*d)**2/pi)/(R*R)
        inverse_distance.append(near+far)
    gradient_main=iv.sqrt(2/pi+3*constant**2+6*pi*pi*sum(inverse_distance))
    gradient_square_cutoff=iv.sqrt(2*tail)
    gradient_hole_cutoff=iv.sqrt(sum(4*pi*K(d,2*rho) for d in norms))
    gradient_total=gradient_main+gradient_square_cutoff+gradient_hole_cutoff
    rough_heat_error=cutoff_difference+pi*iv.sqrt(eps)*gradient_total
    a=pi*eps;vacuum_error=iv.sqrt(2*a*a/((1+a)*(1+2*a)))
    total=2*vacuum_error+rough_heat_error
    def pack(v):
        return dict(display=str(v),dyadic_endpoints=[list(x) for x in v._mpi_])
    values=dict(tail=tail,mass_holes=mass_holes,cutoff_difference=cutoff_difference,
        remainder=remainder,gradient_constant=constant,gradient_main=gradient_main,
        gradient_total=gradient_total,rough_heat_error=rough_heat_error,
        vacuum_error=vacuum_error,total_matrix_element_cutoff_bound=total)
    return dict(dps=dps,epsilon=epsilon,radius=radius,values={k:pack(v) for k,v in values.items()},
        below_point0011=bool(total<iv.mpf('.0011')),
        assumptions=dict(h=1,inner_square=8,outer_square=9,local_radius='.5',radial_split='.6',
            distinct_poles=poles,gradient_poles=gradient_poles),
        scope='Analytic upper bound, interval evaluation with mpmath.iv; library not formally verified. No D_epsilon quadrature certification.')

def main():
    root=Path(__file__).resolve().parents[1];base=root/'results/stage4/heat-matrix'
    base.mkdir(exist_ok=True);out=base/'cutoff-bound.json'
    if out.exists():raise RuntimeError('Do not overwrite')
    rows=[evaluate(p,e,'0.0001') for p in [40,70] for e in ['0.00000001','0.0000000001','0.000000000001']]
    data=dict(rows=rows,versions=dict(python=platform.python_version(),mpmath=m.__version__),
        hashes={'code/heat_cutoff_bound.py':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    out.write_text(json.dumps(data,indent=2))
    for r in rows:print(r['dps'],r['epsilon'],r['values']['total_matrix_element_cutoff_bound']['display'],flush=True)
if __name__=='__main__':main()
