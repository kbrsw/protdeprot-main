from scipy.optimize import minimize

def function(a):
    x,y,z = a
    eq = abs(370*8+((150*x+90*y+110*z)*8)-8900)
    return eq

boundx = (0,3)
boundy = (0,3)
boundz = (0,3)

# opt = minimize(function, (0,0), method = 'SLSQP', bounds=(boundx, boundy))
opt = minimize(function, (0,0,0), bounds=(boundx, boundy, boundz))
print("x = ", opt.x[0])
print("y = ", opt.x[1])
print("z = ", opt.x[2])