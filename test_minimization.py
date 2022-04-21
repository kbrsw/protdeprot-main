from scipy.optimize import minimize

def function(z):
    x,y = z
    eq = ((150*x+137*y)*28-5000)**2
    return eq

boundx = (0,10)
boundy = (0,10)
# opt = minimize(function, (0,0), method = 'SLSQP', bounds=(boundx, boundy))
opt = minimize(function, (0,0), bounds=(boundx, boundy))
print("x = ", opt.x[0])
print("y = ", opt.x[1])