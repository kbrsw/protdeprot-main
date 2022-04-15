import sympy as sp
import math
import mpmath
from scipy.optimize import root, newton
a = 1
#x = sp.symbols('x')
def func1(x):
    f = 2 * x*a - 7
    return f
solution = newton(func1, 10, fprime=lambda x: 2)
diff = 3-solution


i = 0
while diff > 0.0000001:
    def func(x):
        f = 2*(x*a-i/100)-6
        return f
    solution = newton(func, 10, fprime = lambda x: 2)
    i+=1
    diff = 3-solution
    print(solution)