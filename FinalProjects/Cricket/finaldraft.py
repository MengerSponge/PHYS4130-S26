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
from FinalProject_Header import *

# initial conditions
damping = 0.2
vertical_offset = 0.3
t_span = (0,20)
initial_conditions, solutions = [], []

# loop over multiple magnet configurations
magnet_counts = [3, 4, 5, 6, 7, 8]

for idx, N in enumerate(magnet_counts):
    print("")
    ax = axes[idx]

    magnets = generate_polygon_magnets(N)

    restoring = 0.0
    #restoring = 0.6 if N <= 5 else 0.3
    magnetic_strength = 1.2 if N <= 5 else 2.0

    # plot the trajectories for N magnets
    plot_random_trajectories(magnets, f"Trajectories ({N} Magnets)")
    plt.show()

# generating fractals
print("")
print("Generating fractals...")
fig, axes = plt.subplots(2, 3, figsize=(15, 10))  # 6 plots total
axes = axes.flatten()
for idx, N in enumerate(magnet_counts):

    ax = axes[idx] 
    magnets = generate_polygon_magnets(N)
    restoring = 0.0
    magnetic_strength = 1.2 if N <= 5 else 2.0

    # generate fractal
    basins, shading = generate_fast_fractal(magnets, res=400)

    # color pallette
    base_colors = [
        '#e6194b', '#3cb44b', '#ffe119', '#4363d8',
        '#f58231', '#911eb4', '#46f0f0', '#f032e6',
        '#bcf60c', '#fabebe', '#008080', '#e6beff',
        '#9a6324', '#fffac8', '#800000', '#aaffc3',
        '#808000', '#ffd8b1', '#000075', '#808080'
    ]

    colors = base_colors[:N]  # take first N colors
    custom_cmap = ListedColormap(colors)

    # plot
    ax.imshow(basins, cmap=custom_cmap,
              extent=[-1.5, 1.5, -1.5, 1.5],
              origin='lower')

    ax.imshow(shading, cmap='bone',
              extent=[-1.5, 1.5, -1.5, 1.5],
              origin='lower', alpha=0.25)

    # plot magnets
    for m in magnets:
        ax.scatter(m[0], m[1], color='white', s=20)

    ax.set_title(f"{N} Magnets")
    ax.set_xticks([])
    ax.set_yticks([])

plt.suptitle("Magnetic Pendulum Fractals (Various Magnet Configurations)", fontsize=16)
plt.show()

#################################################################
# END FINAL PROJECT
#################################################################
