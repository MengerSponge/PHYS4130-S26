'''
    Project name: Numeric Integration Startup
    Author: Emma Krebs
    Final due date: 2/17/26
    Project description: This project has the student use their knowledge of numerical integration methods
'''

import numpy as np
from scipy.special import legendre
import scipy as sp
import pandas as pd
import math
import matplotlib.pyplot as plt

# ------------------ Function Definitions ------------------

'''
    Definition: trapezoid_function
    Parameters: func (Given function), a (Starting point of subdivision), b (Ending point of subdivision),
        and n (the number of subintervals).
    Description: This function loops over a range of n subintervals and continously sums using the trapezoid 
        equation and our given parameters to get the approximate answer of the integral for n subintervals.
'''
def trapezoid_function(func, a, b, n):

    h = (b - a) / n # Width of our subinterval
    sum = 0
    # Loop for number of subintervals
    for i in range(n):
        sum += (1 / 2) * h * (func(a + i * h) + func(a + i * h + h))

    return sum


'''
    Definition: error_calc
    Parameters: current (Current sum) and true (true value)
    Description: Calculates the error percentage between the given sum and the true value of integral.
'''
def error_calc(current, true):
    return round((abs(current - true) / true) * 100, 5)


'''
    Definition: approximator
    Parameters: f (Given function), start (Starting point of subdivision), end (Ending point of subdivision),
        true_value (Correct value of integral that we will use for approximation), and value of sig_fig we
        are interested in.
    Description: **********
'''
def approximator(f, start, end, true_value, sig_fig):

    n = 1 # Starting number of subintervals
    sums = [] # Empty array to keep track of the sum of the function using n integrals
    parameter = False # Parameter used to know when to end loop once we have approximated close to the answer

    while(parameter == False):

        current_sum = trapezoid_function(f, start, end, n)
        error = error_calc(current_sum, true_value)
        sums.append((n, current_sum, str(error) + '%'))
        
        if error <= sig_fig:
            parameter = True
        else:
            n *= 2

    return sums


'''
    Definition: given_function
    Parameters: x (When called it will plug in the x value)
    Description: Simplifies writing out the function when passed to the trapezoid function
'''
def given_function(x):
    return math.sin(math.sqrt(100*x)) ** 2


'''
    Definition: u (gaussian_quadrature)
    Parameters: x (Value it needs to convert), a (start), and b (end)
    Description: Maps a, b to [-1, 1] interval. By doing so, you can convert an x in that range.
'''
def u(x, a, b):
    return (2*x - a - b) / (b - a)


'''
    Definition: du (gaussian_quadrature)
    Parameters: a (start), and b (end)
    Description: Replaces du so any integral can use dx as the respective variable.
'''
def du(a, b):
    return 2 / (b - a)


'''
    Definition: Legendre_multiplier
    Parameters: 
    Description: 
'''
def Legendre_multiplier(i, j, x_values):
    p1 = legendre(i+1)
    p2 = legendre (j+1)

    u_x_values = u(x_values, min(x_values), max(x_values))

    y = p1(u_x_values) * p2(u_x_values) # Evaulated function for plotting
    integrand = p1 * p2 # Unevaulated function for integrating

    scale = du(min(u_x_values), max(u_x_values)) 

    integrated_value, _ = sp.integrate.quad(integrand, min(u_x_values), max(u_x_values))
    integrated_value = integrated_value * scale 

    return y, integrated_value

# ------------------ Main Body ------------------

# For trapezoid section
results = approximator(given_function, 0, 2, 1.00570254283, 0.00005)
# Create table
print()
pd.set_option('display.precision', 10)
df = pd.DataFrame(results, columns=['Subintervals', 'Summation', 'Error(%)'])
print(df)

# For Gaussian Quadrature (Test)
print("\nTest for Gaussian Quadrature with a = 1, b = 7, and x = 5")
print("The starting point 1 becomes: ", u(1, 1, 7))
print("The ending point 7 becomes: ", u(7, 1, 7))
print("The test variable 5 becomes: ", u(5, 1, 7))
print() # Cheap way to add space

# For Legendre Polynomials
x = np.linspace(-1, 1, 200) 
fig, axes = plt.subplots(4, 4, figsize=(10,10)) # Our subplot grid

for i in range(4):
    y_1 = legendre(i+1)(x)

    for j in range(4):
        y_2 = legendre(j+1)(x)
        axes[i, j].plot(x, y_1, label=f'P{i+1}')
        axes[i, j].plot(x, y_2, label=f'P{j+1}')
        title = f"P{i+1}, P{j+1}, P{i+1}*P{j+1}"
        axes[i, j].set_title(title, size=8)
        y_3, value = Legendre_multiplier(i, j, x)
        axes[i, j].plot(x, y_3, label=f'P{i+1}*P{j+1}')

        axes[i, j].set_xlabel('X')
        axes[i, j].set_ylabel('Legendre Value')
        axes[i, j].legend(fontsize=5)

        print(f"P{i+1}*P{j+1} is {value} and simplified {round(value, 5)}")

plt.tight_layout()
plt.show()

roots_array = []
weights_array = []
values_array = []
for N in range(4):
    roots, weights = sp.special.roots_legendre(N+1)
    roots_array.append(roots)
    weights_array.append(weights)
    values_array.append(N+1)

dg = pd.DataFrame({'P(x)': values_array,'Roots': roots_array, 'Weights': weights_array})
print(dg)

# ------------------ Extension ------------------

# This section is for extension question 1


