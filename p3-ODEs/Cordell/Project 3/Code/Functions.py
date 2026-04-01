import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import RK45
from scipy.integrate import solve_ivp
from scipy.integrate import odeint

# variables
k = 1.0
m = 1.0
b = 0.2
tau =1.0


# functions
def decay_deriv(N,t):
    return -(1./tau) * N


def RK2_Solver(N0, tmin, tmax, nts, decay_deriv):
    N_array = np.zeros(nts)                                                     # array to hold number of nuclei
    t_array = np.linspace(tmin, tmax, nts, endpoint=False)                      # array holds the time points 
    dt = t_array[1] - t_array[0]                                                # dt = time step length  
    N_array[0] = N0                                                             # Initial number of nuclei
    for it in range(0, len(t_array)-1 ):                                        # loop over time steps
        t  = t_array[it]                                                        
        N_h = N_array[it] + (dt/2 * decay_deriv(N_array[it], t))                # sub-step 1 for RK2
        N_array[it+1] = N_array[it] + (dt * decay_deriv(N_h, t + dt/2))         # sub-step 2 for RK2
    return t_array, N_array

def Euler_solver(N_initial, tmin, tmax, nts, deriv):
    N = np.zeros(nts+1)
    t = np.linspace(tmin, tmax, nts+1)
    dt = t[1] - t[0]
    N[0] = N_initial
    for it in range(0,nts):
        N[it+1] = N[it] + dt * deriv(N[it], t[it])
    return t, N

def diffeq_solver_from_scipy(N0, tmin, tmax, nts, decay_deriv):
    t = np.linspace(tmin, tmax, nts, endpoint=False)  
    N = odeint(decay_deriv, N0, t)
    return t, N

def SHO_deriv(xv_array, t):    # x_array is 2D with x, v 
    x, v = xv_array
    dxdt = v
    dvdt = -k/m * x
    return [dxdt, dvdt]

def SHO_solver_RK2(x0, v0, tmin, tmax, nts, SHO_deriv):
    x_array = np.zeros(nts)                                                     # array to hold position
    v_array = np.zeros(nts)                                                     # array to hold velocity
    t_array = np.linspace(tmin, tmax, nts, endpoint=False)                      # array holds the time points 
    dt = t_array[1] - t_array[0]                                                # dt = time step length  
    x_array[0] = x0                                                             # Initial position
    v_array[0] = v0                                                             # Initial velocity
    for it in range(0, len(t_array)-1 ):                                        # loop over time steps
        t  = t_array[it]                                                        
        x_h = x_array[it] + (dt/2 * SHO_deriv([x_array[it], v_array[it]], t)[0])                # sub-step 1 for RK2
        v_h = v_array[it] + (dt/2 * SHO_deriv([x_array[it], v_array[it]], t)[1])                # sub-step 1 for RK2
        x_array[it+1] = x_array[it] + (dt * SHO_deriv([x_h, v_h], t + dt/2)[0])         # sub-step 2 for RK2
        v_array[it+1] = v_array[it] + (dt * SHO_deriv([x_h, v_h], t + dt/2)[1])         # sub-step 2 for RK2
    return t_array, x_array, v_array

def SHO_solver_ODEINT(x0, v0, tmin, tmax, nts, SHO_deriv):
    t = np.linspace(tmin, tmax, nts, endpoint=False)  
    x_v_array = odeint(SHO_deriv, [x0, v0], t)
    x_array = x_v_array[:, 0]
    v_array = x_v_array[:, 1]
    return t, x_array, v_array

def SHO_solver_Euler(x0, v0, tmin, tmax, nts, SHO_deriv):
    x_array = np.zeros(nts)                           
    v_array = np.zeros(nts)                           
    t_array = np.linspace(tmin, tmax, nts, endpoint=False)
    dt = t_array[1] - t_array[0]                        
    x_array[0] = x0                                     
    v_array[0] = v0                                     
    
    for it in range(0, nts-1):
        x_array[it+1] = x_array[it] + dt * SHO_deriv([x_array[it], v_array[it]], t_array[it])[0]
        v_array[it+1] = v_array[it] + dt * SHO_deriv([x_array[it], v_array[it]], t_array[it])[1]
    
    return t_array, x_array, v_array


def SHO_deriv_damped(x_array, t):
    x, v = x_array
    dxdt = v
    dvdt = -k/m * x - b/m * v
    return [dxdt, dvdt]


def H(x, v):
    T = 0.5 * m * v**2
    V = 0.5 * k * x**2
    return T + V

def fun(t, y):
    x, v = y
    dxdt = v
    dvdt = -k/m * x - b/m * v
    return [dxdt, dvdt]

def fun_SHO(t, y):
    x, v = y
    dxdt = v
    dvdt = -k/m * x
    return [dxdt, dvdt]

def A_verlet_SHO(x_array, v_array):
        return -k/m * x_array                       # acceleration (x'') for SHO

def A_verlet_damped(x_array, v_array):
        return -k/m * x_array - b/m * v_array       # acceleration (x'') for damped SHO

def verlet_solver(x0, v0, tmin, tmax, nts, deriv):
    x_array = np.zeros(nts)                                                     # array to hold position
    v_array = np.zeros(nts)                                                     # array to hold velocity                                               
    t_array = np.linspace(tmin, tmax, nts, endpoint=False)                      # array holds the time points 
    dt = t_array[1] - t_array[0]                                                # dt = time step length  
    x_array[0] = x0                                                             # Initial position
    v_array[0] = v0                                                             # Initial velocity
    for it in range(0, len(t_array)-1 ):                                        # loop over time steps
        # Algorithm for Verlet method 
        x1 = x_array[it] + v_array[it]*dt + 0.5 * deriv(x_array[it], v_array[it]) * dt**2
        v1 = v_array[it] + 0.5 * (deriv(x_array[it], v_array[it]) + deriv(x1, v_array[it])) * dt
        x_array[it+1] = x1
        v_array[it+1] = v1


    return t_array, x_array, v_array

def analytical_SHO(x0, v0,t):
    omega = np.sqrt(k/m)
    A = np.sqrt(x0**2 + (v0/omega)**2)
    phi = np.arctan(v0/(omega*x0))
    return A * np.cos(omega*t - phi)

def analytical_damped_SHO(x0, v0,t):
    omega0 = np.sqrt(k/m)
    gamma = b/(2*m)
    omega_d = np.sqrt(omega0**2 - gamma**2)
    A = np.sqrt(x0**2 + ((v0 + gamma*x0)/omega_d)**2)
    phi = np.arctan((v0 + gamma*x0)/(omega_d*x0))
    return A * np.exp(-gamma*t) * np.cos(omega_d*t - phi)

def relative_error(numerical, analytical):
    return np.abs(numerical - analytical) / np.abs(analytical)

def A_verlet_SHO(x_array, v_array):
        return -k/m * x_array                       # acceleration (x'') for SHO

def A_verlet_damped(x_array, v_array):
        return -k/m * x_array - b/m * v_array       # acceleration (x'') for damped SHO