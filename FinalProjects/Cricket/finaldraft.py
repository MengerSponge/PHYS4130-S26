'''
Filename: finaldraft.ipynb
Written by: Cricket Bergner
Date: 04/27/26
'''

# import libraries
import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import solve_ivp as intg
import random as ra
from matplotlib.colors import ListedColormap # for the fractal mapping
from matplotlib import animation
from FinalProject_Header import *

# initial conditions
damping = 0.2
vertical_offset = 0.3
t_span = (0,50)
initial_conditions, solutions = [], []

# loop over multiple magnet configurations
magnet_counts = [3, 4, 5, 6, 7, 8]

print("")
print("Generating motion plots...")

for idx, N in enumerate(magnet_counts):
    magnets = generate_polygon_magnets(N)

    restoring = 0.0
    #restoring = 0.6 if N <= 5 else 0.3
    magnetic_strength = 1.2 if N <= 5 else 2.0

    # plot the trajectories for N magnets
    plot_random_trajectories(magnets, f"Trajectories ({N} Magnets)")
    plt.savefig(f"trajectories_{N}_magnets.png", bbox_inches='tight', dpi=300)
    plt.show()

# generating fractals
print("")
print("Generating fractals...")
fig, axes = plt.subplots(2, 3, figsize=(15, 10))  # 6 plots total
axes_nn = axes.flatten()
for idx, N in enumerate(magnet_counts):

    ax = axes_nn[idx] 
    magnets = generate_polygon_magnets(N)
    restoring = 0.0
    magnetic_strength = 1.2 if N <= 5 else 2.0

    # generate fractal
    basins, shading = generate_fast_fractal(magnets, res=400)

    # color pallette
    base_colors = [
        '#582c83', '#FFD100', '#0fc7f2', "#49f866",
        '#fb0cce', '#3e5c1f', '#3e5cf6', '#fe5c1f',
        '#fc0000', '#fabebe', '#008080', '#e6beff',
    ]

    colors = base_colors[:N]  # take first N colors
    custom_cmap = ListedColormap(colors)

    # plot
    ax.imshow(basins, cmap=custom_cmap,
              extent=[-1.5, 1.5, -1.5, 1.5],
              origin='lower',
              vmin=0, vmax=N-1, 
              interpolation='nearest',
              aspect='equal') 

    ax.imshow(shading, cmap='bone',
              extent=[-1.5, 1.5, -1.5, 1.5],
              origin='lower', alpha=0.25,
              aspect='equal')

    # plot magnets
    for m in magnets:
        ax.scatter(m[0], m[1], color='white', s=20)

    ax.set_title(f"{N} Magnets")
    ax.set_xticks([])
    ax.set_yticks([])

plt.suptitle("Magnetic Pendulum Fractals (Various Magnet Configurations)", fontsize=16)
plt.savefig("magnetic_pendulum_fractals_grid.png", bbox_inches='tight', dpi=300)
plt.show()

# generating animations

print("")
print("Generating animations...")

configs = [3, 4, 5] # number of magnets
step_sizes = [512, 2048] # to show how step sizes effects the accuracy

for N in configs:
    magnets = generate_polygon_magnets(N)
    restoring = 0.05
    magnetic_strength = 1.2

    # initial conditions are the same, but this is set up so you could change them later
    states = [
        [-1.2, 1.20, 0.0, 0.0], # RK45
        [-1.2, 1.20, 0.0, 0.0], # DOP853
        [-1.2, 1.20, 0.0, 0.0]  # Heun
    ]

    for steps in step_sizes:
        filename = f"{N}_magnets_{steps}_steps"
        animate_integrators(magnets, states, steps, filename)

print("")
print("Code finished running!")

###############################################################################################################
# END FINAL PROJECT
###############################################################################################################
