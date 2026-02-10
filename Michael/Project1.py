
import numpy as np
import pylab as py
import scipy as sp
roots, weights = sp.special.roots_legendre(N)

def trapezoid(f,a,b,N):
    h = (b-a)/N #Interval size
    mysum = 0
    
    for i in range(1,N): #should go from 1 to N-1
        mysum = mysum + float(f(a+(i)*h)*h) 
    return mysum +(h/2)*(f(a) + f(b))

def f(x):
    return (np.sin(np.sqrt(100*x)))**2

#Trapezoidal rule stuff
TrueVal = trapezoid(f,0,2,5000000) 

print('\n',f"{'Intervals':<12} {'Approx Value':>16} {'Error' :>18}")
for i in [2**n for n in range(18)]:
    val = trapezoid(f,0,2,i)
    print(f"{i:10d} {val:19.10f} {np.abs(TrueVal - val):18.10f}")
print()

#Gaussan Quadrature stuff
