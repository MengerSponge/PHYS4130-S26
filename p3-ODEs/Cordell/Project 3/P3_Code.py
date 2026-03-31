#############################################################

# going to pull everything from the nb and organize it like this

#############################################################

# 1) Select SHO or SHO with damping
# 2) Select numerical method 
# 3) Plot phase space
# 4) Plot Energy vs Time
# 5) Recurse or end program

#############################################################


from Functions import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.integrate import odeint

print("Select system:")
print("1. Simple Harmonic Oscillator (SHO)")
print("2. Damped Harmonic Oscillator")
system_choice = int(input("Enter your choice (1 or 2): "))

##Fix RK4

print("\nSelect numerical method:")
print("1. Euler's Method")
print("2. RK2")
print("3. RK4")
print("4. Verlet Integration")
print("5. Scipy's ODEINT")
method_choice = int(input("Enter your choice (1, 2, 3, 4, or 5): "))


if system_choice == 1:
    x0 = float(input("Enter initial position: "))
    v0 = float(input("Enter initial velocity: "))
    tmin = 0
    tmax = float(input("Enter final time: "))
    nts = int(input("Enter number of time steps: "))
    if method_choice == 1:
        t, x, v = SHO_solver_Euler(x0, v0, tmin, tmax, nts, SHO_deriv)
    elif method_choice == 2:
        t, x, v = SHO_solver_RK2(x0, v0, tmin, tmax, nts, SHO_deriv)
    elif method_choice == 3:
        t, x, v = solve_ivp(fun, (tmin, tmax), [x0, v0], t_eval=np.linspace(tmin, tmax, nts, endpoint=False), method='RK45').t, solve_ivp(fun, (tmin, tmax), [x0, v0], t_eval=np.linspace(tmin, tmax, nts, endpoint=False), method='RK45').y[0], solve_ivp(fun, (tmin, tmax), [x0, v0], t_eval=np.linspace(tmin, tmax, nts, endpoint=False), method='RK45').y[1]
    elif method_choice == 4:
        t, x, v = verlet_solver(x0, v0, tmin, tmax, nts, A_verlet_SHO)
    elif method_choice == 5:
        t, x, v = SHO_solver_ODEINT(x0, v0, tmin, tmax, nts, SHO_deriv)
        pass
    else:
        print("Invalid method choice for simple harmonic oscillator.")

if system_choice == 2:
    x0 = float(input("Enter initial position: "))
    v0 = float(input("Enter initial velocity: "))
    tmin = 0
    tmax = float(input("Enter final time: "))
    nts = int(input("Enter number of time steps: "))
    if method_choice == 1:
        t, x, v = SHO_solver_Euler(x0, v0, tmin, tmax, nts, SHO_deriv_damped)
    elif method_choice == 2:
        t, x, v = SHO_solver_RK2(x0, v0, tmin, tmax, nts, SHO_deriv_damped)
    elif method_choice == 3:
        t, x, v = solve_ivp(fun, (tmin, tmax), [x0, v0], t_eval=np.linspace(tmin, tmax, nts, endpoint=False), method='RK45').t, solve_ivp(fun, (tmin, tmax), [x0, v0], t_eval=np.linspace(tmin, tmax, nts, endpoint=False), method='RK45').y[0], solve_ivp(fun, (tmin, tmax), [x0, v0], t_eval=np.linspace(tmin, tmax, nts, endpoint=False), method='RK45').y[1]
    elif method_choice == 4:
        t, x, v = verlet_solver(x0, v0, tmin, tmax, nts, A_verlet_damped)
    elif method_choice == 5:
        t, x, v = SHO_solver_ODEINT(x0, v0, tmin, tmax, nts, SHO_deriv_damped)
    else:
        print("Invalid method choice for damped harmonic oscillator.")

print("\nSelect plot to display:")
print("1. Phase Space")
print("2. Energy vs Time")
plot_choice = int(input("Enter your choice (1 or 2): "))

if plot_choice == 1:
    plt.plot(x, v)
    plt.xlabel("Position")
    plt.ylabel("Velocity")
    if system_choice == 1:
        plt.title("Phase Space of Simple Harmonic Oscillator")
    else:
        plt.title("Phase Space of Damped Harmonic Oscillator")
    plt.grid()
    plt.show()

elif plot_choice == 2:
    energy = H(x, v)
    plt.plot(t, energy)
    plt.xlabel("Time")
    plt.ylabel("Energy")
    if system_choice == 1:
        plt.title("Energy vs Time for Simple Harmonic Oscillator")
    else:
        plt.title("Energy vs Time for Damped Harmonic Oscillator")
    plt.grid()
    plt.show()

## relative error


# Target error? nts incrementation ?(not sure if powers of 2 are the best))o