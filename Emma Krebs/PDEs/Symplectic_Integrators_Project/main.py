'''
    Project name: Symplectic Integrations
    Subfolder: main.py
    Author: Emma Krebs
    Final due date: 4/1/26
    Project description: 
'''


import ODE_Methods
import matplotlib.pyplot as plt


# Let us create a dictionary for our possible situations.
# The key is our angular frequency value and our dampening term. The other values of initial 
# conditions, max time, and number of steps will remain the same.

# All of our different situations we can use
situations = {
        'w = 0.5, damp = 0': [1, 1, 150, 0.5, 0, 3000],
        'w = 1.0, damp = 0': [1, 1, 150, 1, 0, 3000],
        'w = 1.5, damp = 0': [1, 1, 150, 1.5, 0, 3000],
        'w = 1.0, damp = 0.25': [1, 1, 150, 1.0, 0.25, 3000],
        'w = 1.0, damp = 0.5': [1, 1, 150, 1.0, 0.5, 3000]
}

# Our different methods of integrations. 
methods = {
    'Verlet Symplectic': ODE_Methods.Verlet_symplectic,
    'Odeint': ODE_Methods.Odeint_solver,
    'RK45': ODE_Methods.RK45_solver,
    'Analytic': ODE_Methods.Harmonic_oscillator
}


# What do the phase spaces look like for these integrators for our different situations.
for method_name, method_func in methods.items():
    plt.figure()
    for value in situations.values():
        x_array, p_array, t_array = method_func(*value)
        found_key = ODE_Methods.Find_key(situations, value)
        plt.plot(x_array, p_array, label=f"{found_key}")
    
    plt.axis('equal')
    plt.title(f'Phase Space of {method_name}')
    plt.xlabel('x Values')
    plt.ylabel('p Values')
    plt.legend()
    plt.grid()

plt.show()


# ----------------- Considering Error -------------------------

energy = [] # Stores all the energy arrays for each method
# How do these graphs compare against each other? Let us take the w = 1 and damp = 0 cases.
for method_name, method_func in methods.items():
    x_array, p_array, t_array = method_func(1, 1, 150, 1, 0, 3000)
    plt.plot(x_array, p_array, label=f"Method: {method_name}")

    kinetic_energy, potential_energy, total_energy = ODE_Methods.Total_energy(x_array, p_array, 1)
    energy.append([kinetic_energy, potential_energy, total_energy, t_array])


found_key = ODE_Methods.Find_key(situations, [1, 1, 150, 1, 0, 3000])
plt.axis('equal')
plt.title(f'Phase Space of N = 3000')
plt.xlabel('x Values')
plt.ylabel('p Values')
plt.legend()
plt.grid()
plt.show()

# Compare Total Energy. Loop through energy and grab each array for a method of energy arrays 
for i in range(len(energy)):
    # Keep track which is which
    if i==0:
        method_name = 'Symplectic'
    elif i==1:
        method_name = 'Odeint'
    elif i==2:
        method_name = 'RK45'
    else:
        method_name = 'Analytic'

    N = energy[i]
    plt.plot(N[3], N[2], label=method_name) # Plots t_array and total_energy array

plt.legend()
plt.title('Total Energy vs Time')
plt.xlabel('Time')
plt.ylabel('Total Energy')
plt.axis([0, 25, 0.9985, 1.001]) # Zoomed in on axis
plt.show()

# Relative Error
relative_error_symp = []
relative_error_ode = []
relative_error_rk = []

analytic_total_array = energy[3][2]
symp_total_array = energy[0][2]
odient_total_array = energy[1][2]
rk_total_array = energy[2][2]

for analytic, symp, odient, rk in zip(analytic_total_array, symp_total_array, odient_total_array, rk_total_array):
    error_symp = ODE_Methods.Relative_error(analytic, symp)
    error_ode = ODE_Methods.Relative_error(analytic, odient)
    error_rk = ODE_Methods.Relative_error(analytic, rk)

    relative_error_symp.append(error_symp)
    relative_error_ode.append(error_ode)
    relative_error_rk.append(error_rk)

# Zoomed out over entire time span
plt.plot(energy[3][3], relative_error_symp, label='Symp')
plt.plot(energy[3][3], relative_error_ode, label='Odeint')
plt.plot(energy[3][3], relative_error_rk, label='RK45')
plt.legend()
plt.axis([0, 25, 0, 0.001]) # Zoomed in to compare first 25 seconds
plt.xlabel('Time')
plt.ylabel('Relative Error')
plt.title('Relative Error for Methods Energy vs Time')
plt.show()

# Compare Time Steps and how that affects the energy drift

N_values = [100, 500, 1000, 3000]

for method_name, method_func in methods.items():

    # We don't care about the analytic for time steps because it shouldn't change
    if method_name == 'Analytic':
        continue 

    plt.figure()

    for N in N_values:

        x_array, p_array, t_array = method_func(1, 1, 150, 1, 0, N)

        kinetic, potential, total = ODE_Methods.Total_energy(x_array, p_array, 1)

        energy_drift = []

        for t in total:
            energy_drift.append(t - total[0])

        h = 150 / N # Time step where 150 is our time duration

        plt.plot(t_array, energy_drift, label=f"N={N}")

    plt.title(f"Energy Drift vs Time for {method_name}")
    plt.xlabel("Time (s)")
    plt.ylabel("Energy Drift (E[t] - E[0])")
    plt.legend()
    plt.show()


# -------------------- Extensions -----------------------------

# def Virial_theorem():


