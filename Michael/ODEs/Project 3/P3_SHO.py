import numpy as np
import matplotlib.pyplot as plt
from P3_Methods import RK45, LSODA, VelVerlet, Analytic, Yoshida

# Notes ----------------------------------------------------------------------------------------------------

#As a convenience, the solvers are configured to solve: dx/dt = f(t, x,y), dy/dt = g(t, x,y)
#Where x and y may be numpy arrays. This works because numpy functions are vectorized.

#For ease of use and reusability of code, all the methods are implimented in the same way:
# t, X, Y = method(X0, Y0, tmin, tmax, nts, du_dt)
# t is the array of times
# X is the array of numerical solutions to x(t)
# Y is the array of numiercal solutions to y(t)
# X0 = X initial condition
# Y0 = Y initial condition
# tmin = starting time
# tmax = finishing time
# nts = numer of points to compute at between tmin and tmax
# du_dt is the vector of derivatives , du_dt = [dx_dt, dy_dt] = [f(t,x,y), g(t,x,y)]



# Derivatives ---------------------------------------------------------------------------------------------------------

# Simple harmonic oscillator. Take x = 0 as the equillibrim point.
# The derivative functions need to return a list [dx_dt, dy_dt]

k = 1.1 # N/m        (spring constant)
m = 1 # kg         (mass)
c = 0.5 # N/(m/s) (damping term strength)

def SHO(t, u): #we need to express the system using u = [x, y]. (we need t as an argument, but its not always used)
    x,y = u #u is the phase space state of the system.

    #Put the derivative for x and y in here
    #x is the position. y is the velocity
    dx_dt = y
    dy_dt = -(k/m)*x
    return [dx_dt, dy_dt]

def SHO_damped(t, u): 
    x,y = u

    #Put the derivative for x and y in here
    #x is the position. y is the velocity
    dx_dt = y
    dy_dt = -(k/m)*x  - (c/m)*y
    return [dx_dt, dy_dt]

def Energy(X,Y): #Energy at different values of time
    return  0.5*m*(Y**2) + 0.5*k*(X**2)



# PLOTTING THE DIFFERENT THINGS *****************************************************************************
# Time Stuff and Initial Conditions -----------------------------------------------------------------------------------------------

tmin = 0 #s start time
tmax = 25 #s end time
nts = 150 #number of points between tmin and tmax

X0 = 1 #m Initial positon
Y0 = 0 #m/s Initial velocity

E0 = 0.5*m*(Y0**2) + 0.5*k*(X0**2)



# UNDAMPED CALCULATIONS -------------------------------------------------------------------------------------
# RK45 Calculations ---------------------------------------------------------------------------------------------------------

solution = RK45([X0], [Y0], tmin, tmax, nts, SHO)
t = solution[0]
XRK_Undamped = solution[1][0]
YRK_Undamped = solution[1][1]

plt.plot(XRK_Undamped,YRK_Undamped, label = "RK4(5)")
plt.title("RK4(5) Phase Space (Undamped)")
plt.xlabel("Position (m)")
plt.ylabel("Velocity (m/s)")
plt.show()

plt.plot(t, Energy(XRK_Undamped, YRK_Undamped))
plt.axhline(y = E0, color = 'r', linestyle = '-')
plt.title("RK4(5) Energy vs Time (Undamped)")
plt.xlabel("Time (s)")
plt.ylabel("Energy (J)")
plt.show()

# LSODA Calculations ---------------------------------------------------------------------------------------------------------

solution = LSODA([X0], [Y0], tmin, tmax, nts, SHO)
t = solution[0]
XLSODA_Undamped = solution[1][0] 
YLSODA_Undamped = solution[1][1]

plt.plot(XLSODA_Undamped,YLSODA_Undamped, label = "LSODA")
plt.title("LSODA Phase Space (Undamped)")
plt.xlabel("Position (m)")
plt.ylabel("Velocity (m/s)")
plt.show()

plt.plot(t, Energy(XLSODA_Undamped, YLSODA_Undamped))
plt.axhline(y = E0, color = 'r', linestyle = '-')
plt.title("LSODA Energy vs Time (Undamped)")
plt.xlabel("Time (s)")
plt.ylabel("Energy (J)")
plt.show()



# Velocity Verlet Calculations ---------------------------------------------------------------------------------------------------------

solution = VelVerlet([X0], [Y0], tmin, tmax, nts, SHO)
t = solution[0]
XVerlet_Undamped = solution[1][0]
YVerlet_Undamped = solution[1][1]

plt.plot(XVerlet_Undamped, YVerlet_Undamped, label = "VelVerlet")
plt.title("Velocity Verlet Phase Space Plot (Undamped)")
plt.xlabel("Position (m)")
plt.ylabel("Velocity (m/s)")
plt.show()

plt.plot(t, Energy(XVerlet_Undamped, YVerlet_Undamped))
plt.axhline(y = E0, color = 'r', linestyle = '-', label = "True Energy")
plt.title("Velocity Verlet Energy vs Time (Undamped)")
plt.xlabel("Time (s)")
plt.ylabel("Energy (J)")
plt.show()



# DAMPED CALCULATIONS  ---------------------------------------------------------------------------------------------------------
# RK45 Calculations ---------------------------------------------------------------------------------------------------------

solution = RK45([X0], [Y0], tmin, tmax, nts, SHO_damped)
t = solution[0]
XRK_Damped = solution[1][0]
YRK_Damped = solution[1][1]

plt.plot(XRK_Damped, YRK_Damped, label = "RK4(5)")
plt.title("RK4(5) Phase Space Plot (Damped)")
plt.xlabel("Position (m)")
plt.ylabel("Velocity (m/s)")
plt.show()

plt.plot(t, Energy(XRK_Damped, YRK_Damped ))
plt.title("RK4(5) Energy vs Time (Damped)")
plt.xlabel("Time (s)")
plt.ylabel("Energy (J)")
plt.show()



# LSODA Calculations ---------------------------------------------------------------------------------------------------------

solution = LSODA([X0], [Y0], tmin, tmax, nts, SHO_damped)
t = solution[0]
XLSODA_Damped = solution[1][0]
YLSODA_Damped = solution[1][1]

plt.plot(XLSODA_Damped, YLSODA_Damped, label = "LSODA")
plt.title("LSODA Phase Space (Damped)")
plt.xlabel("Position (m)")
plt.ylabel("Velocity (m/s)")
plt.show()

plt.plot(t, Energy(XLSODA_Damped, YLSODA_Damped))
plt.title("LSODA Energy vs Time (Damped)")
plt.xlabel("Time (s)")
plt.ylabel("Energy (J)")
plt.show()



# Velocity Verlet Calculations ---------------------------------------------------------------------------------------------------------

solution = VelVerlet([X0], [Y0], tmin, tmax, nts, SHO_damped)
t = solution[0]
XVerlet_Damped = solution[1][0] 
YVerlet_Damped = solution[1][1]

plt.plot(XVerlet_Damped, YVerlet_Damped, label = "VelVerlet")
plt.title("Velocity Verlet Phase Space (Damped)")
plt.xlabel("Position (m)")
plt.ylabel("Velocity (m/s)")
plt.show()

plt.plot(t, Energy(XVerlet_Damped, YVerlet_Damped))
plt.title("Velocity Verlet Energy vs Time (Damped)")
plt.xlabel("Time (s)")
plt.ylabel("Energy (J)")
plt.show()


# END PLOTTING ***********************************************************************************************

# ERROR CALCULATIONS/COMPARISONS ***********************************************************************************************
# This system is analytically solvable for the under damped, critically damped, and over damped cases
# All are handled with the function named "Analytic" which is imported from methods

# Checking a numeric solution with the analytic solution -----------------------------------------------------------

t = np.linspace(tmin, tmax, nts)

plt.plot(t, Analytic(X0, Y0, c, k, m, t), label = "Analytic Solution")
plt.plot(t, XLSODA_Damped, label = "Numeric", linestyle = "--")
plt.legend()
plt.show()



# Relative Error Computations and Plotting --------------------------------------------------------------------------

#redefine these here so that I don't have to keep scrolling to change tmin and tmax

X0 = 4
Y0 = 0

k = 0.75
m = 1.22

tmin = 0
tmax = 30

def RelErr(Analytic, Numeric): #This is cleaner than typing out this formula a bunch imo
    return np.abs((Numeric - Analytic)/Analytic)

nts = range(50, 200)

errRK =[]
errLSODA = []
errVerlet = []
errYoshida = []
dt_array = []

X_Analytic = Analytic(X0, Y0, 0, k, m, tmax)
for i, n in enumerate(nts):

    solution = RK45([X0], [Y0], tmin, tmax, n, SHO)
    X = solution[1][0]
    errRK.append(RelErr(X_Analytic, X[len(X) - 1]))

    solution = LSODA([X0], [Y0], tmin, tmax, n, SHO)
    X = solution[1][0]
    errLSODA.append(RelErr(X_Analytic, X[len(X) - 1]))

    solution = VelVerlet([X0], [Y0], tmin, tmax, n, SHO)
    X = solution[1][0]
    errVerlet.append(RelErr(X_Analytic, X[len(X) - 1]))

    solution = Yoshida([X0], [Y0], tmin, tmax, n, SHO)
    X = solution[1][0]
    errYoshida.append(RelErr(X_Analytic, X[len(X) - 1]))

    t = solution[0]
    dt_array.append(t[1] - t[0])

plt.plot(dt_array, errRK, label = "RK4(5)")
plt.plot(dt_array, errLSODA, label = "LSODA")
plt.plot(dt_array, errVerlet, label = "Velocity Verlet")
plt.plot(dt_array, errYoshida, label = "Yoshida")

plt.xlabel("Step Size (s)")
plt.ylabel("Relative Error")
plt.title("Relative Error vs Step Size")
plt.legend()
plt.show()
