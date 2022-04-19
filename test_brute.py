from scipy.optimize import brute
import numpy as np


def f1(z):
    x,y = z
    return np.sin(2 / x + y)


rranges = (slice(1, 10, 1), slice(1, 10, 1))
optParA = brute(f1, rranges, full_output=True)
print(" x = ", optParA[0], "y = ", optParA[1])
