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
# conditions, step size, and number of steps will remain the same.
situations = {
        'w = 0.5, damp = 0': [1, 1, 150, 0.5, 0, 300],
        'w = 1.0, damp = 0': [1, 1, 150, 1, 0, 300],
        'w = 1.5, damp = 0': [1, 1, 150, 1.5, 0, 300],
        'w = 1.0, damp = 0.25': [1, 1, 150, 1.0, 0.25, 300],
        'w = 1.0, damp = 0.5': [1, 1, 150, 1.0, 0.5, 300]
}

methods = {
    # 'Verlet Symplectic': ODE_Methods.Verlet_symplectic,
    'Odeint': ODE_Methods.Odeint_solver,
    'RK45': ODE_Methods.RK45_solver,
    'Analytic': ODE_Methods.Harmonic_oscillator
}

for method_name, method_func in methods.items():
    for value in situations.values():
        x_array, p_array, t_array = method_func(*value)
        found_key = ODE_Methods.Find_key(situations, value)
        plt.plot(x_array, p_array, label=f"Method name: {method_name} for {found_key}")
    
plt.axis('equal')
plt.title('Phase Space')
plt.xlabel('x Values')
plt.ylabel('p Values')
plt.legend()
plt.grid()
plt.show()
