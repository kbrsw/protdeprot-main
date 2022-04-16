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

optxy = fsolve(optfunc,(xin,yin))
print(optxy)

residue1 = optxy[0]-optxy[1]-5
print("new residue = ", residue1)

residue = residue1 # !!! give another name that we use in cycling

#----- found new A and B that agree optimized X and Y

def fA(A):
    eq1A = A*optxy[0]-5
    return eq1A

def fB(y):
    eq2B = B*optxy[1]+6
    return eq2B

An = fsolve(fA, A)
Bn = fsolve(fB, B)

print("new A = ", An)
print("new B = ", Bn)

A = An
B = Bn
