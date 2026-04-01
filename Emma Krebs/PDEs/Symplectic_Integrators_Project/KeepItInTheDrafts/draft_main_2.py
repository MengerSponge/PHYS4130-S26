'''
    Project name: Symplectic Integrations
    Subfolder: main.py
    Author: Emma Krebs
    Final due date: 2/26/26
    Project description: 
'''


import ODE_Methods
import matplotlib.pyplot as plt


# ------------------ Symplectic Euler Solver ------------------

# Let us create a dictionary for our possible situations.
# The key is our angular frequency value and our dampening term. The other values of initial 
# conditions, step size, and number of steps will remain the same.
situations_symplectic = {
        'w = 0.5, damp = 0': [1, 1, 0.5, 0.5, 300, 0],
        'w = 1.0, damp = 0': [1, 1, 1, 0.5, 300, 0],
        'w = 1.5, damp = 0': [1, 1, 1.5, 0.5, 300, 0],
        'w = 1.0, damp = 0.25': [1, 1, 1, 0.5, 300, 0.25],
        'w = 1.0, damp = 0.5': [1, 1, 1, 0.5, 300, 0.5]
}

# Loop through the parameters in our dictionary and plot them
for value in situations_symplectic.values():
    x_array, p_array, t_array = ODE_Methods.Symplectic_Euler(*value)
    found_key = ODE_Methods.find_key(situations_symplectic, value)
    plt.plot(x_array, p_array, label=found_key)

# Graph information
plt.axis('equal')
plt.title('Simple Harmonic Oscillator Momentum vs Position (Symplectic Euler)')
plt.xlabel('X Values')
plt.ylabel('P Values')
plt.legend()
plt.grid()
plt.show()


# -------------------- RK45 Solver ----------------------------

situations_RK45 = {
    'w = 0.5, damp = 0': [(1, 1), 0.5, 0, 0, 150, 0.5],
    'w = 1.0, damp = 0': [(1, 1), 1, 0, 0, 150, 0.5],
    'w = 1.5, damp = 0': [(1, 1), 1.5, 0, 0, 150, 0.5],
    'w = 1.0, damp = 0.25': [(1, 1), 1, 0.25, 0, 150, 0.5],
    'w = 1.0, damp = 0.5': [(1, 1), 1, 0.5, 0, 150, 0.5]
}

for value in situations_RK45.values():
    x_array, p_array, t_array = ODE_Methods.RK45_solver(*value)
    found_key = ODE_Methods.find_key(situations_RK45, value)
    plt.plot(x_array, p_array, label=found_key)

# Graph information
plt.axis('equal')
plt.title('Simple Harmonic Oscillator Momentum vs Position (Symplectic Euler)')
plt.xlabel('X Values')
plt.ylabel('P Values')
plt.legend()
plt.grid()
plt.show()

# ODE_Methods.RK45_solver()

# ------------------- Odeint Solver ---------------------------

situations_Odeint = {
    'w = 0.5, damp = 0': [(1, 1), 0.5, 0, 0, 150, 75],
    'w = 1.0, damp = 0': [(1, 1), 1, 0, 0, 150, 75],
    'w = 1.5, damp = 0': [(1, 1), 1.5, 0, 0, 150, 75],
    'w = 1.0, damp = 0.25': [(1, 1), 1, 0.25, 0, 150, 75],
    'w = 1.0, damp = 0.5': [(1, 1), 1, 0.5, 0, 150, 75]

}

for value in situations_Odeint.values():
    t_array, N_array = ODE_Methods.Odeint_solver(*value)
    found_key = ODE_Methods.find_key(situations_Odeint, value)
    plt.plot(N_array[:,0], N_array[:,1], label=found_key)

# Graph information
plt.axis('equal')
plt.title('Simple Harmonic Oscillator Momentum vs Position (Symplectic Euler)')
plt.xlabel('X Values')
plt.ylabel('P Values')
plt.legend()
plt.grid()
plt.show()

# ----------------- Considering Error -------------------------

# How do these graphs compare against each other? Let us take the w = 1 and damp = 0 cases.
t_array_ode, N_array_ode = ODE_Methods.Odeint_solver((1, 1), 1, 0, 0, 150, 75) # Odeint
x_array_rk, p_array_rk, t_array_rk = ODE_Methods.RK45_solver((1, 1), 1, 0, 0, 150, 0.5) # RK45
x_array_symp, p_array_symp, t_array_symp = ODE_Methods.Symplectic_Euler(1, 1, 1, 0.5, 300, 0) # Symplectic 

plt.plot(N_array_ode[:,0], N_array_ode[:,1], label='Odeint')
plt.plot(x_array_rk, p_array_rk, label='RK45')
plt.plot(x_array_symp, p_array_symp, label='Symplectic')

plt.title('Comparing Three Methods using w=1 and damp=0')
plt.legend()
plt.grid()
plt.show()

# Compare total mechanical energy
kinetic_symp, potential_symp, total_symp = ODE_Methods.Total_energy(x_array_symp, p_array_symp, 1)
kinetic_ode, potential_ode, total_ode = ODE_Methods.Total_energy(N_array_ode[:,0], N_array_ode[:,1], 1)
kinetic_rk, potential_rk, total_rk = ODE_Methods.Total_energy(x_array_rk, p_array_rk, 1)

# Get analytic solution
x_array_analytic, p_array_analytic, t_array_analytic = ODE_Methods.Harmonic_oscillator(150, 75, (1,1), 1, 0)
kinetic_analytic, potetnial_analytic, total_analytic = ODE_Methods.Total_energy(x_array_analytic, p_array_analytic, 1)

# COmpare
plt.plot(t_array_ode, total_ode, label='Odeint')
plt.plot(t_array_rk, total_rk, label='RK')
plt.plot(t_array_symp, total_symp, label='Symplectic')
plt.plot(t_array_analytic, total_analytic, label='Analytic')
plt.legend()
plt.show()

# Relative Error
relative_error_symp = []
relative_error_ode = []
relative_error_rk = []
for analytic, estimate_symp, estimate_ode, estimate_rk in zip(total_analytic, total_symp, total_ode, total_rk):
    error_symp = ODE_Methods.Relative_error(analytic, estimate_symp)
    error_ode = ODE_Methods.Relative_error(analytic, estimate_ode)
    error_rk = ODE_Methods.Relative_error(analytic, estimate_rk)

    relative_error_symp.append(error_symp)
    relative_error_ode.append(error_ode)
    relative_error_rk.append(error_rk)

plt.plot(t_array_analytic, relative_error_symp)
# plt.plot(t_array_analytic, relative_error_ode)
# plt.plot(t_array_analytic, relative_error_rk)
plt.show()

# -------------------- Extensions -----------------------------


