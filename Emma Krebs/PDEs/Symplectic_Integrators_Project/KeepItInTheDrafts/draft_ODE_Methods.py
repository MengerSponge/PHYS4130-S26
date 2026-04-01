'''
    Project name: Symplectic Integrations
    Subfolder: ODE Methods
    Author: Emma Krebs
    Final due date: 2/26/26
    File description: This file is to be used in the main.py. It includes the main three ODE's needed for 
                    comparison of the harmonic oscillator. This inlcudes a syplectic integration, RK45, 
                    and Odeint solvers.
'''


import numpy as np
import math
from scipy.integrate import RK45, odeint


def Harmonic_oscillator(t, h, y, w, damp):
    t_array = np.linspace(0, t, h)
    x_0 = y[0]
    p_0 = y[0]
    x_array = []
    p_array = []

    for t in t_array:
        x = x_0 * math.cos(w*t) + p_0 * math.sin(w*t) / w
        p = p_0 * math.cos(w*t) - w*x_0*math.sin(w*t)
        x_array.append(x)
        p_array.append(p)
    
    return x_array, p_array, t_array

        
def Relative_error(analytic, estimated):
    return np.abs(analytic-estimated) / analytic


def Harmonic_deriv(t, y, w, damp):
    
    """
    Hammonic deriviation function for the other ODEs to use. Incorporates variables x and p
    as free values and returns how they change at that point in time.

    Args:
        t (int/float): Time (not used in this definition but passed to use in future ODEs).
        y (tuple): Starting position. Starting momentum. What the system looks like.
        w (int/float): Angular frequency of harmonic oscillator.
        damp (int/float): Damping term on oscillator (if there is one)

    Returns:
        tuple: Functions for how dxdt and dpdt change in time.
    """
    x, p = y

    dxdt = p
    dpdt = -(w**2)*x - damp*p

    return [dxdt, dpdt]


def Symplectic_Euler(x_0, p_0, w, h, steps, damp):

    """
    Modification of the Euler method for solving Hamiltonains. Yields better results than Euler
    by similarily using time steps but essentially flipped the order of the velocity update
    and the position update, computing the new velocity of the particle before we compute it's new position.
    In this function we use a symplectic euler method such that we calculate the new momentum before 
    the new position.

    Args:
        x_0 (int/float): Initial x position.
        p_0 (int/float): Initial momentum position.
        w (int/float): Angular frequency for harmonic oscillator.
        h (int/float): Size of step (how much are we changing over)
        steps (int): Number of iterations for function.
        damp (int/float): Damping term (if there is one).

    Returns:
        Three arrays: x_array and p_array of harmonic oscillator for a certain time/steps. Also include t_array from steps.
    """

    x = x_0
    p = p_0

    x_array = [x_0]
    p_array = [p_0]
    t_array = []
    h_total = 0 # To count the current time step
    t_array.append(h_total)

    # Iterations, update momentum first then position
    for i in range(steps):

        p = p - h*((w**2)*x + damp*p)
        x = x + h*p

        x_array.append(x)
        p_array.append(p)
        h_total += h # Increase total time by the step
        t_array.append(h_total)
    
    return x_array, p_array, t_array


def RK45_solver(y_0, w, damp, tmin, tmax, step):

    """
    Calls the RK45 solver from scipy.integrate which is the Runge–Kutta–Fehlberg method, an 
    algorithm numerical analysis for the numerical solution of ODE's. 

    Args:
        y_0 (tuple): Initial conditions.
        w (int/float): Angular frequency of harmonic oscillator.
        damp (int/float): Dampening term.
        tmin (int/float): Minimum time value.
        tmin (int/float): Maximum time value.
        step (int/float): Value of steps using for the time function
        
    Returns:
        Three arrays: Arrays for x values and p values for time. t_array for time is returned as well. 
    """

    # Object that is used to generate arrays
    solver = RK45(lambda t, y_0: Harmonic_deriv(t, y_0, w, damp), tmin, y_0, tmax, max_step=step)

    x_array = []
    p_array = []
    t_array = []

    # Repeat until solver stops running (when it hits the tmax after n steps)
    while solver.status == 'running':
        x_array.append(solver.y[0]) # Current value for x
        p_array.append(solver.y[1]) # Current value for p
        t_array.append(solver.t) # Current time
        solver.step()

    return x_array, p_array, t_array


def Odeint_solver(y_0, w, damp, tmin, tmax, number_steps):

    """
    Calls the odeint function from scipy.integrate. odeint solves a system of ordinary 
    differential equations using lsoda from the FORTRAN library odepack.

    Args:
        y_0 (tuple): Initial conditions.
        w (int/float): Angular frequency of harmonic oscillator.
        damp (int/float): Dampening term.
        tmin (int/float): Minimum time value.
        tmin (int/float): Maximum time value.
        number_steps (int/float): Number of steps

    Returns:
        tuple: Returns t as the t_array and N as an tuple of arrays, which incorporates 
        all the x and p values. 
    """

    t = np.linspace(tmin, tmax, number_steps)
    N = odeint(Harmonic_deriv, y_0, t, args=(w, damp), tfirst=True)
    
    return t, N


def find_key(dict, value):
    """
    Grabs the key for a given value. This is used for labeling of graphs.

    Args: 
        dict (dictionary): Dictionary for the inputs to our solvers.
        value (array of int/floats): The key's value we are searching for.
    
    Returns:
        string: Returns either the key, or if there is no key, None. 
    
    """

    for key, val in dict.items():
        if val == value:
            return key
        
    return None


def Total_energy(array_x, array_p, w):

    """
    This function definition returns the kinetic, potential, and total energy for a given 
    array of x and p (along with their angular frequency) for a harmonic oscillator.

    Args:
        array_x (array): Array of x positions
        array_p (array): Array of p momentums
        w (int.float): Angular frequency

    Returns:
        Three arrays: Kinetic, potential, and total energy arrays
        
    """

    kinetic_energy = []
    potential_energy = []
    total_energy = []
    
    # Assuming m = 1, we can find the kinetic, potential, and total energy through the following methods
    for x, p in zip(array_x, array_p):
        KE_value = (1/2) * p**2
        PE_value = (1/2) * w**2 * x**2 # Since w is angular frequency and m = 1, so k = w**2
        TE_value = KE_value + PE_value

        kinetic_energy.append(KE_value)
        potential_energy.append(PE_value)
        total_energy.append(TE_value)
    
    return kinetic_energy, potential_energy, total_energy
