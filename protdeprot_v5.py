import mpmath
import numpy as np
from sympy import *
from scipy.optimize import root, newton, basinhopping, minimize, fsolve
import sympy as sym
import math as mt
import scipy

sigmaM = float(input("Enter measured charge mineral charge, mol/m2 "))
pH = float(input("Enter pH of solution "))
M = float(input("Enter equilibrium concentration of metal in solution, mol/kg "))
An = float(input("Enter concentration of anion in solution, mol/kg "))
Ns = float(input("Enter a number of surface hydroxyl groups, mol/m2 "))
h = 10 ** (-pH)
oh = (1e-14) / h
n = Ns / 2

# Initial assumptions for constants
# Aain = 10 ** (1.1)
# Abin = 10 ** (3.8)
# Kain = 1e-4
# Kbin = 1.99e-10
Aain = 11
Abin = 11
Kain = 1e-6
Kbin = 1e-6


# -------find som and soh2x at initial Aa, Ka, Ab, Kb
def eq1(x):
    return ((x * h) / ((n - x) * M)) * np.exp(Aain * (x / n)) - Kain


solution1 = fsolve(eq1, 0)
print("som initial =", solution1)


def eq2(y):
    return ((y * oh) / ((n - y) * An)) * np.exp((Abin * y) / (n)) - Kbin


solution2 = fsolve(eq2, 0)
print("soh2x initial = ", solution2)

# ---calculate sigmaCalc
sigmaCalc = solution2 - solution1
residue = sigmaCalc - sigmaM
print("residue = ", residue)



i = 0
pp = []
while abs(residue) > 1e-12:
    # ------- Optimizing som and soh2x taking into account condition sigmaCalc = sigmaM
    def optfunc(x):

        op1 = (solution2[0] - solution1[0] - x[1] + x[0])**2
        op2 = (x[1] - x[0] - sigmaM)**2
        return (op1+op2)

    bound1 = (1e-12,1)
    bound2 = (1e-12,1)
    optConc = minimize(optfunc, (solution1[0], solution2[0]), method = 'Powell', bounds=(bound1, bound2))
    # print(optConc)

    residue1 = optConc.x[1] - optConc.x[0] - sigmaM
    print("new residue = ", residue1)
    pp.append(residue1)
    u = np.array(pp)
    residue = residue1


    # ------Estimate new Aa, Ka, Ab, Kb that agree with opimized som and soh2x

    def f1(z):
        eqAA = (((optConc.x[0] * h) / ((n - optConc.x[0]) * M)) * np.exp(z[0] * (optConc.x[0] / n)) - z[1])**2
        # x[0]=Aa , x[1]=Ka
        return eqAA


    ### !!!!! --- Add boundaries
    boundAa = (10, 300000)
    boundKa = (1e-27, 1e-3)
    optParA = minimize(f1, (Aain, Kain), method='Powell', bounds=(boundAa, boundKa))
    print("new Aa = ", optParA.x[0], "new Ka = ", optParA.x[1])
    Aain = optParA.x[0]
    Kain = optParA.x[1]


    def f2(y):
        eqAB = (((optConc.x[1] * oh) / ((n - optConc.x[1]) * An)) * np.exp((y[0] * optConc.x[1]) / (n)) - y[1])**2
        # y[0] = Ab, y[1]=Kb
        return eqAB


    boundAb = (10, 300000)
    boundKb = (1e-27, 1e-3)
    optParB = minimize(f2, (Abin, Kbin), method='Powell', bounds=(boundAb, boundKb))
    print("new Ab = ", optParB.x[0], "new Kb = ", optParB.x[1])

    Abin = optParB.x[0]
    Kbin = optParB.x[1]

    print ("residue = ", residue)


    solution1[0]=optConc.x[0]
    solution2[0]=optConc.x[1]

    i += 1
