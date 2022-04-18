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
Aain = 0
Abin = 0
Kain = 0
Kbin = 0


def eq1(x):
    return ((x * h) / ((n - x) * M)) * np.exp(Aain * (x / n)) - Kain


solution1 = fsolve(eq1, 0)
print("Equation for Ka and Aa")
print("som =", solution1)


def eq2(y):
    return ((y * oh) / ((n - y) * An)) * np.exp((Abin * y) / (n)) - Kbin


solution2 = fsolve(eq2, 0)
print("Equation for Kb and Ab")
print("soh2x = ", solution2)

# ---- Optimization of som and soh to obtain sigmaCalc = sigmaMeas


# ---Initial guesses are solutions of equations for som and soh2x
a = solution1[0]
b = solution2[0]

x0 = np.array([a, b])

somBound = (1e-09, Ns / 2)
soh2xBound = (1e-09, Ns / 2)


def residue(z):
    f = abs(z[1] - z[0] - sigmaM)
    return f


opt = minimize(residue, x0=x0, method="SLSQP", bounds=(somBound, soh2xBound))
print("som optimized = ", opt.x[0])
print("soh2x optimized = ", opt.x[1])

# ------Substitute optimized som and soh2x in equations for constants and optimizing Ka, Aa, Kb, Ab

som = opt.x[0]
soh2x = opt.x[1]

# ---Equation 1

# --boundaries Aa and Ka
boundAa = (10**1, 10**1.6)
boundKa = (10**(-10.6), 10**(-3.98))

x0 = np.array([Aain, Kain])  # Initial guess


def eqopt1(x):
    return abs(((som * h) / ((n - som) * M)) * np.exp(x[0] * (som / n)) - x[1])


optEq1 = minimize(eqopt1, x0=x0, method="SLSQP", bounds=(boundAa, boundKa))
print("Optimized Aa = ", optEq1.x[0], "Optimized Ka = ", optEq1.x[1])

# ---Equation 2

# --boundaries Ab and Kb
boundAb = (1, 10**3.9)
boundKb = (10**(-27), 10**(-2))

y0 = np.array([Abin, Kbin])  # Initial guess


def eqopt2(y):
    return abs(((soh2x * oh) / ((n - soh2x) * An)) * np.exp((y[0] * soh2x) / (n)) - y[1])


optEq2 = minimize(eqopt2, x0=y0, method="SLSQP", bounds=(boundAb, boundKb))
print("Optimized Ab = ", optEq2.x[0], "Optimized Kb = ", optEq2.x[1])

# ----------------------Take optimized Ka, Aa, Kb, Ab and solve equations again

i = 0
o = []
while abs(opt.x[1] - opt.x[0] - sigmaM) > 1e-12:

    Aain = optEq1.x[0]
    Kain = optEq1.x[1]
    GuessSOM = som


    def eq1(x):
        return ((x * h) / ((n - x) * M)) * np.exp(Aain * (x / n)) - Kain


    solution1 = fsolve(eq1, GuessSOM)
    # print("Equation for Ka and Aa")
    # print("som =", solution1)

    Abin = optEq2.x[0]
    Kbin = optEq2.x[1]
    GuessSOH2X = soh2x


    def eq2(y):
        return ((y * oh) / ((n - y) * An)) * np.exp((Abin * y) / (n)) - Kbin


    solution2 = fsolve(eq2, GuessSOH2X)
    # print("Equation for Kb and Ab")
    # print("soh2x = ", solution2)


    # ----Optimizing SOM and SOH2 to obtain SigmaCalc = SigmaM

    a = solution1[0]
    b = solution2[0]
    x0 = np.array([a, b])

    somBound = (1e-09, Ns / 2)
    soh2xBound = (1e-09, Ns / 2)


    def residue(z):
        f = abs(z[1] - z[0] - sigmaM)
        return f


    opt = minimize(residue, x0=x0, bounds=(somBound, soh2xBound))
    # print("som optimized = ", opt.x[0])
    # print("soh2x optimized = ", opt.x[1])


    # -----Opimization route for parameters 2
    som = opt.x[0]
    soh2x = opt.x[1]
    # ---Equation 1
    x = np.array([0] * 2)

    # --boundaries Aa and Ka
    boundAa = (10 ** 1, 10 ** 1.6)
    boundKa = (10 ** (-10.6), 10 ** (-3.98))

    x0 = np.array([Aain, Kain])  # Initial guess


    def eqopt1(x):
        return abs(((som * h) / ((n - som) * M)) * np.exp(x[0] * (som / n)) - x[1])
    optEq1 = minimize(eqopt1, x0=x0, method="SLSQP", bounds=(boundAa, boundKa))


    # ---Equation 2

    # --boundaries Ab and Kb
    boundAb = (1, 10 ** 3.9)
    boundKb = (10 ** (-27), 10 ** (-2))

    y0 = np.array([Abin, Kbin])  # Initial guess


    def eqopt2(y):
        return abs(((soh2x * oh) / ((n - soh2x) * An)) * np.exp((y[0] * soh2x) / (n)) - y[1])
    optEq2 = minimize(eqopt2, x0=y0, method="SLSQP", bounds=(boundAb, boundKb))

    print("residue = ", solution2[0] - solution1[0] - sigmaM)

    res = soh2x - som - sigmaM
    o.append(res)
    array = np.array(o)

    if abs(opt.x[1] - opt.x[0] - sigmaM) < 1e-12:
        print("Aa = ", optEq1.x[0], "Ka = ", optEq1.x[1], "Ab = ", optEq2.x[0], "Kb = ", optEq2.x[1])

    # if array[i] == array[i+1]:
    #     print("Aa = ", optEq1.x[0], "Ka = ", optEq1.x[1], "Ab = ", optEq2.x[0], "Kb = ", optEq2.x[1])
    #     break

    i += 1
