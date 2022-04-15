import mpmath
import numpy as np
from sympy import *
from scipy.optimize import root, newton, basinhopping, minimize
import time
import math as mt
import scipy


sigmaM = float(input("Enter measured charge mineral charge, mol/m2 "))
pH = float(input("Enter pH of solution "))
M = float(input("Enter equilibrium concentration of metal in solution, mol/kg "))
An = float(input("Enter concentration of anion in solution, mol/kg "))
Ns = float(input("Enter a number of surface hydroxyl groups, mol/m2 "))
h = 10**(-pH)
oh = (1e-14)/h
n = Ns/2

# Initial assumptions for constants
Aain = 12.58
Abin = 398
Kain = 1e-6
Kbin = 1.99e-15

step = 0.5

def eq1(x):
    som = x[0]
    F = np.empty((1))
    F[0]=((x[0]*h)/((n-x[0])*M))*mpmath.exp(Aain*(x[0]/n))
    F[1]=(((sigmaM+x[0])*oh)/((n-sigmaM-x[0])*An))*mpmath.exp(Abin*(sigmaM+x[0])/(n))
    return F
xGuess = np.array([0,0])
solution = newton(eq1, xGuess, maxiter = 5000)
sigmacalc = solution[1]-solution[0]

diff = sigmacalc - sigmaM


print("sigmaCalc = ", sigmacalc)
print("sigmameas=", sigmaM)
print("sigmaCalc-sigmaM = ",  diff)



i = 0
d = []
while abs(diff) > 0.000000001:
    def eqiter(x):
        F = np.empty((2))
        F[0] = ((x[0] * h) / ((n - x[0]) * M)) * mpmath.exp(((Aain) * x[0]) / n) - (Kain-i/5e+5)
        F[1] = ((x[1] * oh) / ((n - x[1]) * An)) * mpmath.exp(((Abin) * x[1]) / n) - (Kbin-i/5e+14)
        return F

    xGuess = np.array([0, 0])
    solution = newton(eqiter, xGuess, maxiter=5000)
    # solution = root(eqiter, xGuess)
    sigmacalc = solution[1] - solution[0]
    diff = sigmacalc - sigmaM
    d.append(diff)
    print("sigmaCalc = ", sigmacalc, "sigmames = ", sigmaM, "error=", diff, "Ka = ", Kain-i/5e+5, "Kb = ", Kbin+i/5e+14)
    u = np.array(d)
# esli umenshenie constant ne podhodit
    if u[i] > u[i-1]:
        i = 0
        print("Another route...")
        time.sleep(2)
        while abs(diff) > 0.000000001:
            def eqiter(x):
                F = np.empty((2))
                F[0] = ((x[0] * h) / ((n - x[0]) * M)) * mpmath.exp(((Aain) * x[0]) / n) - (Kain + i/5e+5)
                F[1] = ((x[1] * oh) / ((n - x[1]) * An)) * mpmath.exp(((Abin) * x[1]) / n) - (Kbin + i/5e+14)
                return F
            xGuess = np.array([0, 0])
            solution = newton(eqiter, xGuess, maxiter=5000)
        # solution = root(eqiter, xGuess)
            sigmacalc = solution[1] - solution[0]
            diff = sigmacalc - sigmaM
            print("sigmaCalc = ", sigmacalc, "sigmames = ", sigmaM, "error=", diff, "Ka = ", Kain + i/5e+5, "Kb = ", Kbin + i/5e+14)
            i += 1

    i += 1


#Output result
print("som concentration ", solution[0])
print("sohx concentrations ",  solution[1])
print("logAa = ", mpmath.log10(Aain))
print("logKa = ", mpmath.log10(Kain-i/5e+5))
print("logKb = ", mpmath.log10(Kbin-i/5e+14))
print("Calculated Ns = ", solution[0]+solution[1], "mol/kg")


