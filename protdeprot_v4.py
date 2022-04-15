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

#---- Optimization of som and soh to obtain sigmaCalc = sigmaMeas
z = np.array([0]*2)
z[0] = solution1[0]
z[1] = solution2[0]


x0 = np.array([0]*2)

somBound = (1e-16,1)
soh2xBound = (1e-16,1)

def residue(z):
    f = abs(z[1] -z[0] - sigmaM)
    return f

opt = minimize(residue, x0=x0, method='SLSQP', bounds=(somBound, soh2xBound))
print("som optimized = ", opt.x[0])
print("soh2x optimized = ", opt.x[1])

#------Substitute optimized som and soh2x in equations for constants and optimizing Ka, Aa, Kb, Ab
som = opt.x[0]
soh2x = opt.x[1]
#---Equation 1
x = np.array([0]*2)
x[0] = Aain
x[1] = Kain

#--boundaries Aa and Ka
boundAa = (0,300)
boundKa = (1e-10, 1e-1)

x0 = np.array([Aain, Kain]) # Initial guess

def eqopt1(x):
    return abs(((som * h) / ((n - som) * M)) * np.exp(x[0] * (som / n)) - x[1])
optEq1 = minimize(eqopt1, x0=x0, method='SLSQP', bounds=(boundAa, boundKa))
print("Optimized Aa = ", optEq1.x[0], "Optimized Ka = ", optEq1.x[1])

#---Equation 2
y = np.array([0]*2)
y[0] = Abin
y[1] = Kbin

#--boundaries Ab and Kb
boundAb = (0,300)
boundKb = (1e-20, 1e-5)

y0 = np.array([Abin, Kbin]) #Initial guess

def eqopt2(y):
    return ((soh2x * oh) / ((n - soh2x) * An)) * np.exp((y[0] * soh2x) / (n)) - y[1]
optEq2 = minimize(eqopt2, x0=y0, method='SLSQP', bounds=(boundAb, boundKb))
print("Optimized Ab = ", optEq2.x[0], "Optimized Kb = ", optEq2.x[1])