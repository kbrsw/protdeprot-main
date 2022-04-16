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
h = 10**(-pH)
oh = (1e-14)/h
n = Ns/2

# Initial assumptions for constants
Aain = 10**(1.1)
Abin = 10**(3.8)
Kain = 1e-4
Kbin = 1.99e-10


def eq1(x):
     return ((x*h)/((n-x)*M))*np.exp(Aain*(x/n)) - Kain
solution1 = fsolve(eq1, 0)
print("Equation for Ka and Aa")
print("som =", solution1)

def eq2(y):
    return ((y * oh) / ((n - y) * An)) * np.exp((Abin * y) / (n)) - Kbin
solution2 = fsolve(eq2, 0)
print("Equation for Kb and Ab")
print("soh2x = ", solution2)

x = np.array([0]*2)
x[0] = solution1[0]
x[1] = solution2[0]

x0 = np.array([0]*2)

somBound = (1e-16,1)
soh2xBound = (1e-16,1)

def residue(x):
    f = abs(x[1] -x[0] - sigmaM)
    return f

# x0 = np.array([x[0], x[1])

print("Residue ", residue(x))

opt = minimize(residue, x0=x0, method='SLSQP', bounds=(somBound, soh2xBound))
print("optimized som = ",opt.x[0])
print("optimized soh2x = ",opt.x[1])
print("optimized residue = ", opt.x[1]-opt.x[0]-sigmaM)


som = opt.x[0]
soh2x = opt.x[1]
#-----
def min(x):
    # x[0] = Aa
    # x[1]=Ka
    # x[2] = Ab
    # x[3] = Kb
    eqGen = abs(((som*h)/((n-som)*M))*np.exp(x[0]*(som/n)) - x[1] - ((soh2x * oh) / ((n - soh2x) * An)) * np.exp((x[2] * soh2x) / (n)) + x[3])
    return eqGen

x0 = np.array([10**(1.35), 1e-4, 10**(3.8), 1.99e-10])

x1bond = (5,100)
x2bond = (1e-13, 1e-6)
x3bond = (10,10000)
x4bond = (1e-18, 1e-6)
opt1 = minimize(min, x0=x0, method='SLSQP', bounds=(x1bond, x2bond, x3bond, x4bond))
print("Aa = ", opt1.x[0],"Ka = ", opt1.x[1], "Ab = ", opt1.x[2], "Kb = ", opt1.x[3])

i=0
while (opt.x[1]-opt.x[0]-sigmaM) > 1e-12:

    # x = opt.x[0]
    # y = opt.x[1]

    def eq1(x):
        return ((x * h) / ((n - x) * M)) * np.exp(opt1.x[0] * (x / n)) - opt1.x[1]
    solution1 = fsolve(eq1, 0)


    def eq2(y):
        return ((y * oh) / ((n - y) * An)) * np.exp((opt1.x[2] * y) / (n)) - opt1.x[3]
    solution2 = fsolve(eq2, 0)

    solution1[0] = som
    solution2[0] = soh2x

    som = x[0]
    soh2x = x[1]

    def residue(x):
        f = abs(x[1] - x[0] - sigmaM)
        return f

        

    print("Residue ", residue(x))

    somBound = (1e-16, 1)
    soh2xBound = (1e-16, 1)

    opt = minimize(residue, x0=x0, method='SLSQP', bounds=(somBound, soh2xBound))
    print("optimized residue = ", opt.x[1] - opt.x[0] - sigmaM)

    # -----
    def min(x):
        # x[0] = Aa
        # x[1]=Ka
        # x[2] = Ab
        # x[3] = Kb
        eqGen = abs(((som * h) / ((n - som) * M)) * np.exp(x[0] * (som / n)) - x[1] - (
                    (soh2x * oh) / ((n - soh2x) * An)) * np.exp((x[2] * soh2x) / (n)) + x[3])
        return eqGen


    x0 = np.array([opt1.x[0], opt1.x[1], opt1.x[2], opt1.x[3]])

    x1bond = (5, 100)
    x2bond = (1e-13, 1e-6)
    x3bond = (10, 10000)
    x4bond = (1e-18, 1e-6)
    opt1 = minimize(min, x0=x0, method='SLSQP', bounds=(x1bond, x2bond, x3bond, x4bond))
    print("Aa = ", opt1.x[0], "Ka = ", opt1.x[1], "Ab = ", opt1.x[2], "Kb = ", opt1.x[3])

    i+=1

# Ka = minimize(residue, [0], method= 'SLSQP', args=(Kain))
# print("optimized Ka", Ka.x[0])
#
# Kb = minimize(residue, [0], method= 'SLSQP', args=(Kbin))
# print("optimized Kb", Kb.x[0])
#
# Aa = minimize(residue, [0], method= 'SLSQP', args=(Aain))
# print("optimized Aa", Aa.x[0])
#
# Ab = minimize(residue, [0], method= 'SLSQP', args=(Abin))
#
# print("optimized Ab", Ab.x[0])
#
# o = []
#
# while soh2x-soh2x-sigmaM > 1e-8:
#
#     v = np.array(o)
#
#     Ka = minimize(residue, [0], method= 'Powell', args=(Ka.x[0]))
#
#
#     Kb = minimize(residue, [0], method= 'Powell', args=(Kb.x[0]))
#
#
#     Aa = minimize(residue, [0], method= 'Powell', args=(Aa.x[0]))
#
#
#     Ab = minimize(residue, [0], method= 'Powell', args=(Ab.x[0]))
#
#
#     def eq1(x):
#         return ((x*h)/((n-x)*M))*np.exp(Aa.x[0]*(x/n)) - Ka.x[0]
#     solution1 = fsolve(eq1, 0)
#
#     def eq2(y):
#         return ((y * oh) / ((n - y) * An)) * np.exp((Ab.x[0] * y) / (n)) - Kb.x[0]
#     solution2 = fsolve(eq2, 0)
#
#
#     som = solution1[0]
#     soh2x=solution2[0]
#
#     sigmaCalc = soh2x-som
#
#     o.append(sigmaCalc)
#
#     print("new residue ", sigmaCalc-sigmaM)
#
#     i+=1

