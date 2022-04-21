import mpmath
import numpy as np
from sympy import *
from scipy.optimize import root, newton, basinhopping, minimize, fsolve, brute, least_squares
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

# Initial guess for constants
Aain = 20
Abin = 3500
Kain = 1e-6
Kbin = 1e-15



residue = 1
x0 = (11,1e-6, 11, 1e-6)
pp1 = []
pp2 = []
pp3 = []
pp4 = []

# -------find som and soh2x at initial Aa, Ka, Ab, Kb
def eq1(x):
    return ((x * h) / ((n - x) * M)) * np.exp(Aain * (x / n)) - Kain


# ----Derivation
x = Symbol('x')
f1 = ((x * h) / ((n - x) * M)) * exp(Aain * (x / n)) - Kain

differ1 = f1.diff(x)
# --- Convert from sympy to numpy expression
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
i = 0
while abs(residue) > 1e-12:

    # ---calculate sigmaCalc
    sigmaCalc = solution2 - solution1
    residue = sigmaCalc - sigmaM
    print("residue = ", residue[0])

    ss1 = solution1[0]
    ss2 = solution2[0]

    # ------- Optimizing som and soh2x taking into account condition sigmaCalc = sigmaM
    def optfunc(z):
        # x[0] = x
        # x[1] = y
        x,y = z
        op1 = abs(ss2 - ss1 - y + x)
        op2 = (y - x - sigmaM)
        return (op1,op2)

    bound1 = (1e-12,1)
    bound2 = (1e-12,1)
    optConc = least_squares(optfunc, (ss1, ss2), bounds = ((1e-16, 1e-16),(1,1)))


    residue = optConc.x[1] - optConc.x[0] - sigmaM

    print("new residue = ", residue)
    print("optimized som ", optConc.x[0], "optimized soh2x ", optConc.x[1])

    # ------Estimate new Aa, Ka, Ab, Kb that agree with opimized som and soh2x

    def f1(z):
        eqAA = abs((((optConc.x[0] * h) / ((n - optConc.x[0]) * M)) * np.exp(z[0] * (optConc.x[0] / n)) - z[1])+(((optConc.x[1] * oh) / ((n - optConc.x[1]) * An)) * np.exp((z[2] * optConc.x[1]) / (n)) - z[3]))
        # x[0]=Aa , x[1]=Ka
        return eqAA


    ### !!!!! --- Add boundaries
    boundAa = (10, 44)
    boundKa = (2.5e-11, 1.05e-4)
    boundAb = (1, 7944)
    boundKb = (1e-27, 1.26e-6)

    optParA = minimize(f1, x0=x0, method='Powell', bounds=(boundAa, boundKa, boundAb, boundKb))
    # optParA = minimize(f1, x0=x0, method = 'SLSQP', bounds=(boundAa, boundKa, boundAb, boundKb))
    x0 = (optParA.x[0], optParA.x[1], optParA.x[2], optParA.x[3])

    print("new Aa = ", optParA.x[0], "new Ka = ", optParA.x[1])
    print("new Ab = ", optParA.x[2], "new Kb = ", optParA.x[3])

    Aain = optParA.x[0]
    Kain = optParA.x[1]
    Abin = optParA.x[2]
    Kbin = optParA.x[3]

    n = optConc.x[0]
    solution1[0] = optConc.x[0]
    solution2[0] = optConc.x[1]

    pp1.append(Aain)
    pp2.append(Kain)
    pp3.append(Abin)
    pp4.append(Kbin)

    i += 1

    u = np.array(pp1)
    v = np.array(pp2)
    e = np.array(pp3)
    w = np.array(pp4)


