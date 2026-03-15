import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import RK45, odeint


def Symplectic_Euler(start_v, start_u, w, h, steps, damp):

    v_array = []
    u_array = []

    v_array.append(start_v)
    u_array.append(start_u)

    prev_v = start_v
    prev_u = start_u

    for i in range(steps):

        curr_v = prev_v + h*prev_u
        curr_u = prev_u - h*((w**2)*curr_v + damp*prev_u)

        v_array.append(curr_v)
        u_array.append(curr_u)

        prev_v = curr_v
        prev_u = curr_u
    
    return v_array, u_array


'''def Harmonic_deriv(intial, t, w, damp):
    
    x_start, p_start = intial
    dxdt = p_start
    dpdt = -(w**2)*x_start - damp*p_start

    return [dxdt, dpdt]'''


def RK45_solver(N_initial, w, damp, tmin, tmax, time_total):
    t = np.linspace(tmin, tmax, time_total)

    solver = RK45(Harmonic_deriv, N_initial, t, args=(w, damp))

    x_val = []
    y_val = []

    while solver.status == 'running':
        
    return x_val, y_val


# Scipy differential solver using odeint
def Odeint_solver(N_initial, w, damp, tmin, tmax, time_total):
    t = np.linspace(tmin, tmax, time_total)

    N = odeint(Harmonic_deriv, N_initial, t, args=(w, damp))
    
    return t, N


# -- Symplectic Euler -- 
x_array_1, p_array_1 = Symplectic_Euler(1, 1, 0.5, 0.05, 300, 0)
x_array_2, p_array_2 = Symplectic_Euler(1, 1, 1, 0.05, 300, 0)
x_array_3, p_array_3 = Symplectic_Euler(1, 1, 1.5, 0.05, 300, 0)
x_array_4, p_array_4 = Symplectic_Euler(1, 1, 1, 0.05, 300, 0.25)

plt.plot(x_array_1, p_array_1, label='w=0.5')
plt.plot(x_array_2, p_array_2, label='w=1')
plt.plot(x_array_3, p_array_3, label='w=1.5')
plt.plot(x_array_4, p_array_4, label='w=1, Damp=0.25')
plt.axis('equal')
plt.title('Simple Harmonic Oscillator Phase Space (Symplectic Euler)')
plt.xlabel('X Values')
plt.ylabel('P Values')
plt.legend()
plt.grid()
plt.show()

# -- RK45 --
# RK45_integrater = RK45()


# -- Odeint --
total_time = 0.05*300
t, N = Odeint_solver([1,1], 1, 0, 0, 1000, 20000)
x_array_o = N[:,0]
p_array_o = N[:,1]
t, N = Odeint_solver([1,1], 1, 0.25, 0, 1000, 20000)
x_array_o_damp = N[:,0]
p_array_o_damp = N[:,1]
# t, N = RK45_solver([1,1], 1, 0.25, 0, 0.05*300, 300)

plt.plot(x_array_o, p_array_o, label='Odeint undamped')
plt.plot(x_array_o_damp, p_array_o_damp, label='Odeint damped')
plt.axis('equal')
plt.title('Simple Harmonic Oscillator Phase Space (Symplectic Euler)')
plt.xlabel('X Values')
plt.ylabel('P Values')
plt.legend()
plt.grid()
plt.show()

