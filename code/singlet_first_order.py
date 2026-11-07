"""Specialization of Lin--Stanford/Streicher connected four-point channels.

h=4 lambda, u_j=2 Js t_j, d=2 for the derivative singlet.
The returned coefficients multiply h, not lambda. Pairwise-disconnected
pieces use full finite-h two-point functions when comparing to exact data.
"""
import numpy as np


def connected_coefficient(times,d=2):
    a,b,c,e=map(float,times)
    g=lambda x:np.cosh(x)**(-2*d)
    direct=-d*d*g(a-b)*g(c-e)*np.tanh(a-b)*np.tanh(c-e)
    nested=-d*d*g(b-c)*g(e-a)*np.tanh(b-c)*np.tanh(e-a)
    crossed=-d*d*g(a-c)*g(b-e)*np.cosh(a+c-b-e)/(np.cosh(a-c)*np.cosh(b-e))
    return dict(direct=direct,nested=nested,crossed=crossed,total=direct+nested+crossed)


def alternating_commutator(u,d=2):
    u=np.asarray(u)
    return 2*d*d*(np.cosh(2*u)+(4*np.tanh(u)**2-1)*np.cosh(u)**(-4*d))


def adjacent_cumulant(u,d=2):
    u=np.asarray(u)
    return d*d*(2*np.tanh(u)**2-1)*np.cosh(u)**(-4*d)
