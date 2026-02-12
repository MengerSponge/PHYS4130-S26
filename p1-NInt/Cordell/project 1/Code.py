
# # # # # # # # # # # # # # # # # 
# # #  Code for Project 1  # # #
# # # # # # # # # # # # # # # # #

# Libraries 
import numpy as np




# # # Trapezoidal Rule # # # 


def leftpoint(f,a,b,N):
    mysum = 0
    h = f
    x_array = np.linspace(a, b, N)
    w = (b-a)/(N-1)
    A_array = h(x_array)*w
    mysum = np.sum(A_array)
    mysum = mysum - A_array[N-1]      # left

    return mysum

def rightpoint(f,a,b,N):
    mysum = 0
    h = f
    x_array = np.linspace(a, b, N)
    w = (b-a)/(N-1)
    A_array = h(x_array)*w
    mysum = np.sum(A_array)
    mysum = mysum - A_array[0]     # right

    return mysum   

def trapezoid(f, a, b, N):
    return 0.5* (leftpoint(f, a, b, N) + rightpoint(f, a, b, N))

def sin(x):
    return (np.sin(sqrt(x)))**2

def sqrt(x):
    return (100*x)**(1/2)


def subintervals(Num_Method, i):
        result = Num_Method(sin, a, b, i)
        error = np.abs(soln - result)
        return (i, result, error)


#def subintervals(Num_Method, N_array):
#    for N in N_array:
#        result = Num_Method(sin, a, b, N)
#        error = np.abs(soln - result)
#        return (N, result, error)
        
soln = 1.00570254282573  # analytic solution via wolfram alpha
a = 0     # lower limit
b = 2     # upper limit

# Start with one single integration slice and work up from there to two, four, eight, and so forth. For each value of the number of slices N
# : your program should print out the number of slices, its estimate of the integral, and its estimate of the error on the integral.

print("TRAPEZOID RULE: ")
N_array = [2**k for k in range(1, 15)]
for i in N_array:
    print("The number of slices, estimate of the integral, and estimated error respectively are:",subintervals(trapezoid, i))


def midpoint(f,a,b,N):        
    h = f
    mysum = 0
    x_array = np.linspace(a, b, N)
    w = (b-a)/(N-1)
    mid_array = [x + (w/2) for x in x_array]
    mid_array = np.array(mid_array)
    A_array = h(mid_array)*w
    mysum = np.sum(A_array)
    mysum = mysum - A_array[N-1]     # 
    return mysum   

def Simpson(f, a, b, N):
    return (1/3 * trapezoid(f, a, b, N)) + (2/3 * midpoint(f, a, b, N))

# Simpson with limits -1 to 1
# PLanning to do the parameterization before sending to Simpson

def sub(u):         # Parameterization: integral from a,b to -1,1 
    x = ((b-a)*u/2)+(a+b)/2
    dx_over_du = 2/(b-a)
    return (1/dx_over_du)*(sin(x))

print("\n")
print("Parameterized Simpson rule result for 1024 subintervals: ", Simpson(sub, -1, 1, 1024))
print("\n")

# Guassian Quadrature

import scipy as sp

def quad(a, b, N):
    roots, weights = sp.special.roots_legendre(N)
    x = ((b-a)*roots/2)+(a+b)/2
    dx_over_du = 2/(b-a)
    return dx_over_du* np.sum(weights*sin(x))

print("GUASSIAN QUADRATURE: ")
N_array = [2**k for k in range(1, 10)]
for i in N_array:
    print ("N: ", i, "; Guassian quadrature result: ", quad(0, 2, i))

# Legendre polynomial plots

import pylab as py
import scipy as sp
P_i = [sp.special.legendre(i) for i in range(1, 5)]
P_j = [sp.special.legendre(i) for i in range(1, 5)]

x = np.linspace(-1, 1, 100)
py.figure(figsize=(14, 14))
for i in range(4):
    for j in range(4):
        Pi = P_i[i](x)
        Pj = P_j[j](x)
        PiPj = Pi*Pj

        axis= py.subplot(4, 4, i*4 +j +1)
        axis.plot(x, Pi, label=f"P{i+1}")
        axis.plot(x, Pj, label=f"P{j+1}")
        axis.plot(x, PiPj, label=f"P{i+1}P{j+1}")

        axis.set_title(f"P{i+1}, P{j+1}, P{i+1}P{j+1}")
        axis.set_xlabel("x")
        axis.set_ylabel("P(x)")
        axis.legend()

py.tight_layout()
py.show()