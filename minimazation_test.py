from sympy import *
import numpy as np
import mpmath
import math
from scipy.optimize import minimize, basinhopping, leastsq, fsolve, root

#---Initial guess for A and B
A = 2
B = 1

#------find x and y at initial A and B
def f1(x):
    eq1 = A*x-5
    return eq1

def f2(y):
    eq2 = B*y+6
    return eq2

xin = fsolve(f1, 1)
print("xin = ", xin)

yin = fsolve(f2, 1)
print("yin = ", yin)

#-----Compare x-y calculated and measured
diffCalc = xin-yin
diffM = 5
residue = diffCalc-diffM
print("residue = ", residue)

# ---- Optimizing X and Y takin into account conditions and given difference between
def optfunc(x):
    op1 = abs(xin[0]-yin[0]-x[0]+x[1])
    op2 = abs(x[0]-x[1]-5)
    return (op1,op2)

optxy = root(optfunc,(xin[0],yin[0]))
print(optxy)

diff= optxy.x[0]-optxy.x[1]
print(diff)
