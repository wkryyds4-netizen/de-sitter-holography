"""BINT four-arc kernel for real times, extending the saved alternating kernel.

Complex bilinear contractions are retained: no absolute square is taken for
general time order. Same-species mutual weight q^4 is included exactly once.
"""
import numpy as np
from singlet_otoc import Kernel


class GeneralKernel(Kernel):
    def evolution(self,t):
        return (self.u*np.exp(1j*t*self.e)[None,:])@self.u.T

    def two(self,t):
        vac=self.u@(np.exp(-1j*t*self.e)*self.u[0,:])
        return float(np.dot(self.rn,abs(vac)**2))

    def four(self,times):
        a,b,c,d=map(float,times)
        arcs=[d-c,c-b,b-a,a-d]
        mats=[self.evolution(s) for s in arcs]
        def uncrossed(shift):
            j=[(i+shift)%4 for i in range(4)]
            left=self.rn*mats[j[2]][:,0];right=self.rn*mats[j[0]][:,0]
            middle=self.evolution(arcs[j[1]]+arcs[j[3]])
            return np.dot(left,middle@right)
        direct=uncrossed(0);nested=uncrossed(1)
        v1=self.weight*mats[0][:,0][self.indices]
        v4=self.weight*mats[3][:,0][self.indices]
        left=mats[2]@v4;right=mats[1]@v1
        terms=self.rn[:,None]*left*right
        reduced=terms.sum();crossed=self.Q*reduced
        total=direct+nested+crossed
        g12=self.two(a-b);g34=self.two(c-d);g23=self.two(b-c);g41=self.two(d-a)
        g13=self.two(a-c);g24=self.two(b-d)
        disc=[g12*g34,g23*g41,g13*g24]
        return dict(times=list(times),arcs=arcs,
            channels=[float(direct.real),float(nested.real),float(crossed.real)],
            channel_imag=[float(direct.imag),float(nested.imag),float(crossed.imag)],
            pair_products=disc,total=float(total.real),connected=float(total.real-sum(disc)),
            crossed_internal_cancellation=float(np.sum(abs(terms))/max(abs(reduced),1e-300)),
            connected_subtraction_condition=float((abs(total)+sum(disc))/max(abs(total-sum(disc)),1e-300)))
