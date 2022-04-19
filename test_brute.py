from scipy.optimize import brute, differential_evolution, shgo, minimize
import numpy as np


def f1(z):
    x,y = z
    return np.sin(2 / x + y)

xr = (1,1000)
yr = (1,1000)

optParA = minimize(f1, (100,100), bounds=(xr,yr))
print(" x = ", optParA.x[0], "y = ", optParA.x[1])
