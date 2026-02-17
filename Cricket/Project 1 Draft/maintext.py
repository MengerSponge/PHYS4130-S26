#########################################################################
#Author: Cricket Bergner
#Date: 02/10/2026
#########################################################################
print("#########################################################################")
print("")
print("")
print("Adaptive Integration")
print("")
print("")
#########################################################################
print("#########################################################################")

# import libraries
import numpy as np
import math as m
from prettytable import PrettyTable as pt

# initialize variables
counter = itl = e = 0
N = 1
ans = 1.005702542825726 #from https://www.integral-calculator.com
t = pt(["Number of Slices", "Integral Estimate", "Estimate of Error"])

# trapezoid rule
def trap(f, a, b, N):
  dx = (b-a)/N
  s = 0 # sum  
  for i in range(N):
    xi = a + (i*dx)
    xip = a + ((i+1)*dx)
    s += ((f(xi) + f(xip))* (dx / 2))
  return s

# while loop to iterate until conditions of problem are met
while counter < 13:
  N *= 2
  counter += 1
  itl = trap(lambda x: np.sin(np.sqrt(100*x))**2, 0, 2, N) # calculate integral
  e =  np.abs(itl - ans) # calculate error
  t.add_row([N, round(itl, 7), round(e, 7)]) # add to table
  
print(t) # table
print("The correct output to the integral is ", ans, ".")
print("It took 8192 intervals for the trapezoid rule to approximate this with an accuracy of 10^-6.")

#########################################################################
print("#########################################################################")
print("")
print("")
print("Gaussian Quadrature")
print("")
print("")

#########################################################################
print("#########################################################################")

# import libraries
import scipy as sp

# define gaussian quadrature
def gauss_quad(f, a, b, N):
    roots, weights = sp.special.root_legendre(N)

    # convert limits of integration
    inv_u = ((b - a) * roots + a + b) / 2
    inv_du = (b - a) / 2

    return(inv_du * np.sum(weights * f(inv_u)))

ans = gauss_quad(lambda x: np.sin(np.sqrt(100*x))**2, 0, 2, 100)
print("Evaluating the integral gives: ", ans)


#########################################################################
print("#########################################################################")
print("")
print("")
print("Subplots")
print("")
print("")
print("#########################################################################")

# import libraries
import scipy as sp
import matplotlib.pyplot as plt
from scipy.special import legendre as l

# initialize variables
x = np.linspace(-1, 1, 400)

# create figure with subplots
f, a = plt.subplots(4, 4, figsize=(16, 16))
f.suptitle('Legendre Polynomials', fontsize=16)

for i in range(4):  # rows
    for j in range(4):  # columns
        ax = a[i, j]
        
        # load Legendre polynomials and their evaluations at x
        # source: https://python4physics.in/program/python/program.php?menu_id=14&submenu_id=2#gsc.tab=0
        Pi = l(i+1)
        Pj = l(j+1)
        Piv = Pi(x) 
        Pjv = Pj(x)
      
        # make the plots
        ax.plot(x, Piv, label=f'P_{i+1}(x)', color='goldenrod', linewidth=1.75)
        ax.plot(x, Pjv, label=f'P_{j+1}(x)', color='purple', linewidth=1.5)
        ax.plot(x, Piv*Pjv, label='Product', color='#78C3F5', linewidth=2)
        
        #label everything for coherency
        ax.set_xlabel('x', fontsize=8)
        ax.set_ylabel('P_n(x)', fontsize=8)
        ax.set_title(f'P_{i+1}, P_{j+1}, P_{i+1}·P_{j+1}', fontsize=10)
        ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)
        ax.tick_params(labelsize=7)

plt.tight_layout()
plt.show()

#########################################################################
print("#########################################################################")
print("")
print("")
print("Extension 1")
print("")
print("")
print("#########################################################################")

# initialize variables
new_counter = 0
new_ans = 6.033977866125206
new_t = pt(["Number of Slices", "Integral Estimate"])
nN = 1

# define necessary functions for Simpson's Rule 
def mid(f, a, b, N): # midpoint rule
  dx = (b-a)/N
  s = 0
  for i in range(N):
    xi = a + (dx/2) + i*dx
    s += (f(xi) * dx)
  return s

def simp(f, a, b, N): # simpson's rule
  dx = (b-a)/N
  s = 0
  s += ((1/3)*trap(f, a, b, N) + (2/3)*mid(f, a, b, N))
  return s

# How many points do you need in your Gaussian quadrature to achieve 10 digits of precision?
while new_counter < 13:
  nN *= 2
  new_counter += 1
  new_itl = gauss_quad(lambda y: (y**2 / np.sqrt(2 - y)), 0, 2, nN)
  new_t.add_row([nN, round(new_itl, 7)]) # add to table

print(new_t)
print("")
print("Running the function beyond this point increases the wait times by a large amount.")
print("Gaussian quadrature is optimized for polynomial-like functions. However, the function")
print("above is not smooth at the endpoint; it has algebraic convergence. This means for large")
print("N values, the big O will increase rapidly, causing unnecessarily long wait times.")
print("Thus, apply a change of variables to fix the issue and cause the integral to behave more like a polynomial.")
print("")

# Apply the change of variable y = 2sin^2(theta)
print("")
print("")
print("Apply change of variable and calculate accuracy using Simpson's Rule.")
print("")

# How many Simpson's rule point do you need to calculate this to 10 significant figures?

# initialize variables
nnew_t = pt(["Number of Slices", "Integral Estimate"])
nnN = 1
nnew_counter = 0

# How many points do you need to achieve 10 digits of precision?
while nnew_counter < 13:
  nnN *= 2
  nnew_counter += 1
  nnew_itl = simp(lambda y: ((16 / np.sqrt(2)) * np.sin(y)**5 ), 0, np.pi/2, nnN)
  nnew_t.add_row([nnN, round(nnew_itl, 12)]) # add to table

print(nnew_t)
print("")
print("The correct output to the integral is ", new_ans, ".")
print("It took ", nnN, "intervals for Simpson's rule to approximate this with an accuracy of 10^-10.")
print("")

# Now use Gaussian quadrature to calculate the integral. 
# How many points do you need to achieve the same precision? 

# initialize variables
nnnew_t = pt(["Number of Slices", "Integral Estimate"])
nnnN = 1
nnnew_counter = 0

# How many points do you need to achieve 10 digits of precision?
while nnnew_counter < 13:
  nnnN *= 2
  nnnew_counter += 1
  nnnew_itl = gauss_quad(lambda y: ((16 / np.sqrt(2)) * np.sin(y)**5 ), 0, np.pi/2, nnnN)
  nnnew_t.add_row([nnnN, round(nnnew_itl, 12)]) # add to table

print(nnnew_t)
print("")
print("The correct output to the integral is ", new_ans, ".")
print("It took ", nnnN, "intervals for Gaussian quadrature to approximate this with an accuracy of 10^-10.")
print("")

#########################################################################
print("#########################################################################")
print("")
print("")
print("Extension 2")
print("")
print("")
print("#########################################################################")



