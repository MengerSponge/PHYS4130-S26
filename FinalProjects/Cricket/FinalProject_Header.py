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
              and assigns each point to a magnet to build the fractal image. Result is stored as an array.
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

    # determine basins
    distances = []
    for m in magnets:
        dist = np.sqrt(np.sum((position - m)**2, axis=-1))
        distances.append(dist)

    basins = np.argmin(distances, axis=0)

    return basins, settled_at

'''
Function name: stop_event
What it does: Triggers if the velocity gets too slow. 
'''
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

'''
Function name: plot_random_trajectories
What it does: Randomly generates initial position and velocity in the x and y directions. 
              Uses the integration function to plot how the motion performs under the random initial conditions for N magnets. 
'''
def plot_random_trajectories(magnets, title):

    initial_conditions, solutions = [], []

    # generate random initial conditions
    for _ in range(3):
        x0 = ra.uniform(-1.5, 1.5)
        y0 = ra.uniform(-1.5, 1.5)
        vx0 = ra.uniform(-0.1, 0.1)
        vy0 = ra.uniform(-0.1, 0.1)
        initial_conditions.append([x0, y0, vx0, vy0])

    # solve trajectories
    for ic in initial_conditions:
        solutions.append(rk45(ic, magnets))

    # plot
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    for i, sol in enumerate(solutions):
        ax = axes[i]
        x = sol.y[0]
        y = sol.y[1]

        ax.plot(x, y, lw=0.7)

        for m in magnets:
            ax.scatter(m[0], m[1], color='red', s=50, zorder=3)

        ic = initial_conditions[i]
        ax.set_title(
            f"x0={ic[0]:.2f}, y0={ic[1]:.2f}\n"
            f"vx0={ic[2]:.2f}, vy0={ic[3]:.2f}"
        )

        ax.set_aspect('equal')
        ax.grid(alpha=0.3)
        ax.set_xlim(-1.5, 1.5)
        ax.set_ylim(-1.5, 1.5)

    plt.suptitle(title)

'''
Function name: heun_non_iterative
What it does: Implements the non-iterative Heun method. Fixed step size integrator.
'''
def heun_non_iterative(initial_state, magnets, dt, steps):

    state = np.array(initial_state, dtype=float)
    trajectory = [state.copy()]

    for _ in range(steps):
        k1 = np.array(equations(0, state, magnets))
        predictor = state + dt * k1
        k2 = np.array(equations(0, predictor, magnets))

        state = state + (dt/2) * (k1 + k2)
        trajectory.append(state.copy())

        # stopping condition (velocity gets too slow)
        speed = np.sqrt(state[2]**2 + state[3]**2)
        if speed < 0.01:
            break

    return np.array(trajectory)

'''
Function name: dop853
What it does: Uses high-order DOP853 integrator from solve_ivp.
'''
def dop853(initial_state, magnets):
    solution = intg(
        fun=lambda t, y: equations(t, y, magnets),
        t_span=t_span,
        y0=initial_state,
        method='DOP853', 
        rtol=1e-7,
        atol=1e-9,
        events=lambda t, y: stop_event(t, y, magnets)
    )
    return solution

'''
Function name: pad_trajectory
What it does: Makes sure that the function does not stop running when the slower integrators have yet to finish. 
'''
def pad_trajectory(traj, target_len):
        current_len = len(traj)
        if current_len < target_len:
            padding = np.tile(traj[-1], (target_len - current_len, 1))
            return np.vstack((traj, padding))
        return traj

'''
Function name: animate_integrators
What it does: Creates animation comparing RK45, Heun, and DOP853. Prints a table of how their final distances vary.
'''
def animate_integrators(magnets, initial_states, steps, filename_prefix):
    dt = (t_span[1] - t_span[0]) / steps

    # generate trajectories with slightly different start points
    rk = rk45(initial_states[0], magnets)
    dp = dop853(initial_states[1], magnets)
    hn = heun_non_iterative(initial_states[2], magnets, dt, steps)

    rk_xy = np.vstack((rk.y[0], rk.y[1])).T
    dp_xy = np.vstack((dp.y[0], dp.y[1])).T
    hn_xy = hn[:, :2]

    # making sure the animations run for all the integrators, even the slower ones. 
    max_len = max(len(rk_xy), len(dp_xy), len(hn_xy))

    rk_xy = pad_trajectory(rk_xy, max_len)
    dp_xy = pad_trajectory(dp_xy, max_len)
    hn_xy = pad_trajectory(hn_xy, max_len)

    # plot
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.set_xlim(-2.0, 2.0)
    ax.set_ylim(-2.0, 2.0)
    ax.set_aspect('equal')
    ax.set_title(f"{len(magnets)} Magnets ({steps} steps)")

    for i, m in enumerate(magnets):
        ax.scatter(m[0], m[1], color='black', s=50, zorder=3)

    rk_line, = ax.plot([], [], color='#71e669', label='RK45 (y0=1.20)')
    dp_line, = ax.plot([], [], color='blue', label='DOP853 (y0=1.20)')
    hn_line, = ax.plot([], [], color='red', label='Heun (y0=1.20)')
    ax.legend(loc='upper right')

    '''
    Function name: update
    What it does: Updates frames with the functions' different positions. 
    '''
    def update(frame):
        rk_line.set_data(rk_xy[:frame, 0], rk_xy[:frame, 1])
        dp_line.set_data(dp_xy[:frame, 0], dp_xy[:frame, 1])
        hn_line.set_data(hn_xy[:frame, 0], hn_xy[:frame, 1])
        return rk_line, dp_line, hn_line

    ani = animation.FuncAnimation(fig, update, frames=max_len - 1, interval=50, blit=True)
    ani.save(f"{filename_prefix}.gif", writer='pillow')
    plt.close()

    # printing the table of distances
    f = max_len - 1
    d_rk_dp = np.linalg.norm(rk_xy[f] - dp_xy[f])
    d_rk_hn = np.linalg.norm(rk_xy[f] - hn_xy[f])
    d_hn_dp = np.linalg.norm(hn_xy[f] - dp_xy[f])
    
    print("")
    print("-" * 40)
    print(f"--- {filename_prefix} ---")
    print(f"Final distance between RK45 and DOP853: {d_rk_dp:.6f}")
    print(f"Final distance between RK45 and Heun:   {d_rk_hn:.6f}")
    print(f"Final distance between Heun and DOP853: {d_hn_dp:.6f}")
    print("-" * 40)

###############################################################################################################
