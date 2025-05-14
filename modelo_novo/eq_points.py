import numpy as np
import sympy as sp
from sympy.utilities.lambdify import lambdify
import scipy.optimize as opt
from sympy import diff
from sympy import sqrt
import decimal as dc
dc.getcontext().prec = 2

g_name = "G"
i_name = "I"
b_name = "B"
G = sp.Symbol(g_name)
I = sp.Symbol(i_name)
B = sp.Symbol(b_name)

RG      = sp.Symbol('RG')
kG      = sp.Symbol('kG')
muG     = sp.Symbol('muG')
alphaI  = sp.Symbol('alphaI')
muI     = sp.Symbol('muI')
alphaB  = sp.Symbol('alphaB')
kB      = sp.Symbol('kB')
muB     = sp.Symbol('muB')

variables = [G,I,B]

def dGdt():
    return RG -kG*I*G -muG*G

def dIdt():
    return alphaI*B -muI*I

def dBdt():
    return alphaB*G*(1000 - B) - muB*B
    
G_eq = dGdt()
I_eq = dIdt()
B_eq = dBdt()

eqMat = sp.Matrix([ G_eq, I_eq, B_eq ])
Mat = sp.Matrix([ G, I, B ])
jacMat = eqMat.jacobian(Mat)
print('Jacobian %s' % jacMat)
print()
eqp = sp.solve((sp.Eq(G_eq, 0.), sp.Eq(I_eq, 0.), sp.Eq(B_eq, 0.)),
    (G, I, B))

for e in eqp:
    print(sp.simplify(e))
    print()
