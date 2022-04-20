from scipy.optimize import newton
import numpy as np
from sympy import *


sigmaM = float(input("Enter measured charge mineral charge, mol/m2 "))
pH = float(input("Enter pH of solution "))
M = float(input("Enter equilibrium concentration of metal in solution, mol/kg "))
An = float(input("Enter concentration of anion in solution, mol/kg "))
Ns = float(input("Enter a number of surface hydroxyl groups, mol/m2 "))
h = 10 ** (-pH)
oh = (1e-14) / h
n = Ns / 2

Aain = 11
Abin = 11
Kain = 1e-6
Kbin = 1e-6


def eq1(x):
     return ((x * h) / ((n - x) * M)) * np.exp(Aain * (x / n)) - Kain


#----Derivation
x = Symbol('x')
f1=((x * h) / ((n - x) * M)) * exp(Aain * (x / n)) - Kain

differ1 = f1.diff(x)
#--- Convert from sympy to numpy expression
deriveq1 = lambdify([x], [differ1], "numpy")


solution1 = newton(eq1, 0, deriveq1)
print(solution1[0])

def eq2(y):
    return ((y * oh) / ((n - y) * An)) * np.exp((Abin * y) / (n)) - Kbin

y = Symbol('y')
f2 = ((y * oh) / ((n - y) * An)) * exp((Abin * y) / n) - Kbin

differ2 = f2.diff(y)
deriveq2 = lambdify([y], [differ2], "numpy")
solution2 = newton(eq2, 0, deriveq2)

print(solution2[0])