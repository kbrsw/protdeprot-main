import mpmath
import numpy as np
from sympy import *
from scipy.optimize import root, newton, basinhopping, minimize, fsolve
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
Aain = 10**(1.1)
Abin = 10**(3.8)
Kain = 1e-6
Kbin = 1.99e-15

step = 0.5

def eq1(x):
    H = np.empty((1))
    H[0] =((x[0]*h)/((n-x[0])*M))*mpmath.exp(Aain*(x[0]/n)) - Kain
    return H
solution1 = fsolve(eq1, 0)
print("Equation for Ka and Aa")
print("som =", solution1[0])

def eq2(y):
    F = np.empty((1))
    F[0] = ((y[0] * oh) / ((n - y[0]) * An)) * mpmath.exp((Abin * y[0]) / (n)) - Kbin
    return F
solution2 = fsolve(eq2, 0)
print("Equation for Kb and Ab")
print("soh2x = ", solution2[0])

sigmaCalc = solution2[0]-solution1[0]
print("Sigma calculated", sigmaCalc)

diff = sigmaCalc-sigmaM
print("Error of sigma = ", diff)

d = [] #array for diff
r = [] #array for Ka

i=0
while abs(diff) > 1e-10:
    def eq1(x):
        H = np.empty((1))
        H[0] = ((x[0] * h) / ((n - x[0]) * M)) * mpmath.exp((Aain * x[0]) / n) - (Kain+i*Kain)
        return H
    solution1 = fsolve(eq1, 0)

    def eq2(y):
        F = np.empty((1))
        F[0] = ((y[0] * oh) / ((n - y[0]) * An)) * mpmath.exp((Abin * y[0]) / n) - (Kbin)
        return F
    solution2 = fsolve(eq2, 0)

    sigmaCalc = solution2[0]-solution1[0]
    diff = sigmaCalc-sigmaM
    d.append(diff)
    r.append(Kain+i*Kain)

    print("Error", diff, "Kb = ", Kbin, "Ka = ", Kain+i*Kain)
    print("sigmaCalc = ", sigmaCalc)

    u = np.array(d)
    v = np.array(r)
    lastKa = v[i-1]

    if abs(u[i])>abs(u[i-1]):
        i = 0
        print("Next route...")
        time.sleep(2)

        o = []  # array for diff
        w = []  #array for Aa
        # while sigmaM != sigmaCalc:
        while abs(diff) > 1e-9:
            def eq1(x):
                H = np.empty((1))
                H[0] = ((x[0] * h) / ((n - x[0]) * M)) * mpmath.exp(((Aain-0.5*i/Aain) * x[0]) / n) - (lastKa)
                return H
            solution1 = fsolve(eq1, 0)

            def eq2(y):
                F = np.empty((1))
                F[0] = ((y[0] * oh) / ((n - y[0]) * An)) * mpmath.exp((Abin * y[0]) / n) - (Kbin)
                return F
            solution2 = fsolve(eq2, 0)

            sigmaCalc = solution2[0] - solution1[0]
            diff = sigmaCalc - sigmaM
            o.append(diff)
            w.append(Aain-0.5*i/Aain)
            print("Error", diff, "Kb = ", Kbin, "Ka = ", lastKa)
            s = np.array(o)
            qq = np.array(w)
            lastAa = qq[i-1]


            if abs(o[i]) > abs(o[i-1]):
                i=0
                o=[]
                w=[]
                while abs(diff)>5e-9:
                    def eq1(x):
                        H = np.empty((1))
                        H[0] = ((x[0] * h) / ((n - x[0]) * M)) * mpmath.exp(((lastAa) * x[0]) / n) - (lastKa)
                        return H
                    solution1 = fsolve(eq1, 0)

                    def eq2(y):
                        F = np.empty((1))
                        F[0] = ((y[0] * oh) / ((n - y[0]) * An)) * mpmath.exp((Abin * y[0]) / n) - (Kbin-i*i*Kbin)
                        return F
                    solution2 = fsolve(eq2, 0)
                    sigmaCalc = solution2[0] - solution1[0]
                    diff = sigmaCalc - sigmaM
                    o.append(diff)
                    w.append(Kbin+50*i*Kbin)
                    print("Error", diff, "Kb = ", Kbin-i*i*Kbin, "Ka = ", lastKa, "Aa= ", lastAa)
                    s = np.array(o)
                    qq = np.array(w)
                    lastKb = qq[i - 1]

                    if abs(o[i])<1e-9:
                        break
                        print("Kb = ", lastKb)

                    if abs(o[i])-abs(o[i-1]) == 0:
                        break
                        print("Kb = ", lastKb)

            i+=1

        i += 1

    i+=1









    # if abs(s[i])>abs(s[i-1]):
    #     i=0
    #     print("Next route!...")
    #     time.sleep(2)
    #     o = [] #array for diff
    #     w = [] #array for Kb
    #
    #     while abs(diff) > 1e-9:
    #         def eq1(x):
    #             H = np.empty((1))
    #             H[0] = ((x[0] * h) / ((n - x[0]) * M)) * mpmath.exp(((lastAa) * x[0]) / n) - (lastKa)
    #             return H
    #         solution1 = fsolve(eq1, 0)
    #
    #         def eq2(y):
    #             F = np.empty((1))
    #             F[0] = ((y[0] * oh) / ((n - y[0]) * An)) * mpmath.exp((Abin * y[0]) / n) - (Kbin+0.01*i*Kbin)
    #             return F
    #         solution2 = fsolve(eq2, 0)
    #         sigmaCalc = solution2[0] - solution1[0]
    #         diff = sigmaCalc - sigmaM
    #         o.append(diff)
    #         w.append(Kbin+0.01*i*Kbin)
    #         print("Error", diff, "Kb = ", Kbin+0.01*i*Kbin, "Ka = ", lastKa, "Aa= ", lastAa)
    #         s = np.array(o)
    #         qq = np.array(w)
    #         lastKb = qq[i - 1]
    #         i += 1



print("Aa = ", lastAa)
print("Ka = ", lastKa)
print("SigmaCalc = ", sigmaCalc)