
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

N_array = [2**k for k in range(1, 15)]
for i in N_array:
    print("The number of slices, estimate of the integral, and estimated error respectively are:",subintervals(trapezoid, i))


