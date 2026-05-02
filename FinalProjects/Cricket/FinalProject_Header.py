'''
Filename: FinalProject_Header.py
Written by: Cricket Bergner
Date: 05/02/26
'''

# import libraries
import numpy as np
from scipy.integrate import solve_ivp as intg
from matplotlib.colors import ListedColormap # for the fractal mapping

'''
Function name: acceleration
What it does: Computes the total acceleration on the pendulum from both the restoring force, damping, and magnetic attractions. 
              Loops over the magnets to find their force contributions to the acceleration in the x and y directions.
'''
def acceleration(x, y, vx, vy, magnets):
    ax = -restoring * x - damping * vx
    ay = -restoring * y - damping * vy

    for mag in magnets:
        dx = mag[0] - x
        dy = mag[1] - y
        r = np.sqrt(dx**2 + dy**2 + vertical_offset**2)
        ax += magnetic_strength * dx / r**3
        ay += magnetic_strength * dy / r**3

    return ax, ay

'''
Function name: equations
What it does: Converts the second-order motion into a system of first-order differential 
              equations for use with the ODE solver. Returns velocity components and acceleration 
              components so the solver can evolve position and velocity over time.
'''
def equations(t, state, magnets):
    x, y, vx, vy = state
    ax, ay = acceleration(x, y, vx, vy, magnets)
    return [vx, vy, ax, ay]

'''
Function name: rk45
What it does: Numerically solves the system of differential equations using an adaptive Runge–Kutta method. 
              Integrates the motion from the initial state until either the time limit is reached or the 
              velocity gets too low.
'''
def rk45(initial_state, magnets):
    solution = intg(
        fun = lambda t, y: equations(t, y, magnets),
        t_span=t_span,
        y0=initial_state,
        method='RK45',
        rtol=1e-5,
        atol=1e-7,
        events=lambda t, y: stop_event(t, y, magnets)
    )
    return solution
  
'''
Function name: find_closest_magnet
What it does: Computes the distance from a given point to each magnet. Determines which magnet is closest.
              Index of closest magnet is returned to classify a basin of attraction.
'''
def find_closest_magnet(x, y, magnets):
    distances = [np.sqrt((x - m[0])**2 + (y - m[1])**2) for m in magnets]
    return np.argmin(distances)

'''
Function name: generate_fast_fractal
What it does: Simulates the motion of many initial conditions simultaneously using vectorized operations instead 
              of individual ODE solvers. It evolves all points forward in time, tracks how quickly they settle, 
              and assigns each point to a magnet to build the fractal image. Result is stored as a 2D array.
'''
def generate_fast_fractal(magnets, res=400, max_iter=400):

    # vectorized grid instead of nested loops
    x = np.linspace(-1.5, 1.5, res)
    y = np.linspace(-1.5, 1.5, res)
    X, Y = np.meshgrid(x, y)

    position = np.stack([X, Y], axis=-1)
    vel = np.zeros_like(position)
    dt = 0.05
    settled_at = np.zeros((res, res))

    for i in range(max_iter):
        accel = -restoring * position - damping * vel

        for m in magnets:
            diff = m - position
            r = np.sqrt(np.sum(diff**2, axis=-1) + vertical_offset**2)
            accel += magnetic_strength * diff / r[..., np.newaxis]**3

        # fast integration (Euler instead of RK45)
        vel += accel * dt
        position += vel * dt
        speed = np.sqrt(np.sum(vel**2, axis=-1))
        settled_at[speed < 0.02] = i

    # Determine basins
    distances = []
    for m in magnets:
        dist = np.sqrt(np.sum((position - m)**2, axis=-1))
        distances.append(dist)

    basins = np.argmin(distances, axis=0)

    return basins, settled_at

# stop event
def stop_event(t, state, magnets):
    vx = state[2]
    vy = state[3]
    speed = np.sqrt(vx**2 + vy**2)
    return speed - 0.01

stop_event.terminal = True
stop_event.direction = -1

'''
Function name: generate_polygon_magnets
What it does: Generates evenly spaced magnet positions around a circle of radius 1 for any number of magnets.
'''
def generate_polygon_magnets(N, R=1.0):
    magnets = [
        np.array([
            R * np.cos(2*np.pi*k/N),
            R * np.sin(2*np.pi*k/N)
        ])
        for k in range(N)
    ]
    return magnets
