from sympy import *
import numpy as np
import mpmath
import math
from scipy.optimize import minimize, basinhopping, leastsq


def func(x):
    eq = abs(2*x[0]+3*x[1]-4000)
    return eq



x0 = np.array([3]*2)
x1bound=(0,5000)
x2bound=(0,6000)



opt = minimize(func, x0=x0, method='SLSQP', bounds=(x1bound, x2bound))

print(opt)
