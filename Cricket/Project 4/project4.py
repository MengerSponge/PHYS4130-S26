'''
Filename: project4.py
Written by: Cricket Bergner
Date: 04/22/2026
'''
# ####################################################################################
# BEGIN PROJECT 4
# ####################################################################################

# import libraries
import numpy as np
from matplotlib import pyplot as plt
import random as ra
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D

# from here on out I have split everything into sections for visual clarity

# ####################################################################################
# SECTION 1: IMPORTANT FUNCTIONS AND INITIAL VARIABLES
# ####################################################################################

# calculated capacity dimension for different DLA functions
def calculate_capacity_dimension(grid):
    # calculates the box-counting dimension
    pixels = np.argwhere(grid > 0)
    if len(pixels) < 2: return 0
    
    # Determine the maximum possible box size based on grid shape
    max_side = max(grid.shape)
    scales = np.unique(np.floor(np.logspace(0, np.log10(max_side / 2), 15)).astype(int))
    scales = scales[scales > 1]
    ns = []
    
    for s in scales:
        bins = (pixels // s).astype(int)
        unique_bins = np.unique(bins, axis=0)
        ns.append(len(unique_bins))
    
    # linear regression
    coeffs = np.polyfit(np.log(scales), np.log(ns), 1)
    return -coeffs[0]

# for 2D DLA, particle wandering
def move_2D(px, py): 
    direction = ra.randint(0, 3)
    if direction == 0: py += 1
    elif direction == 1: py -= 1
    elif direction == 2: px -= 1
    elif direction == 3: px += 1
    return px, py

# for 2D Triangular DLA, particle wandering
def move_2Dt(px, py):
    direction = ra.randint(0, 5)
    if direction == 0: px += 1
    elif direction == 1: px -= 1
    elif direction == 2: py += 1
    elif direction == 3: py -= 1
    elif direction == 4:
        py += 1
        px += (1 if py % 2 != 0 else -1)
    elif direction == 5:
        py -= 1
        px += (1 if py % 2 != 0 else -1)
    return px, py

# for 3D DLA, particle wandering
def move_3D(px, py, pz):
    d = ra.randint(0, 5)
    if d == 0: px += 1
    elif d == 1: px -= 1
    elif d == 2: py += 1
    elif d == 3: py -= 1
    elif d == 4: pz += 1
    elif d == 5: pz -= 1
    return px, py, pz

# if particle nears another particle, it sticks
def sticking(px, py, grid, stickiness):
    if np.any(grid[px-1:px+2, py-1:py+2] > 0):
        return ra.random() <= stickiness
    return False

# if particle nears another particle, it sticks (DLA 2D triangular)
def sticking_2Dt(px, py, grid, stickiness, n):
    dx = 1 if py % 2 != 0 else -1
    neighbors = [(px+1, py), (px-1, py), (px, py+1), (px, py-1),
                 (px+dx, py+1), (px+dx, py-1)]
    for i, j in neighbors:
        if 0 <= i < n and 0 <= j < n and grid[i, j] > 0:
            return ra.random() <= stickiness
    return False

# saves a fram every 50 particles for the 2D DLA animation
def update_2d(frame):
    global spawn, r_max
    
    # run the simulation 50 times per frame 
    particles_per_frame = 50 
    
    for _ in range(particles_per_frame):
        theta = ra.uniform(0, 2*np.pi)
        px, py = int(center + spawn*np.cos(theta)), int(center + spawn*np.sin(theta))
        
        while True:
            px, py = move_2D(px, py)
            
            # Kill radius check
            dist = np.sqrt((px-center)**2 + (py-center)**2)
            if dist > (spawn * kill_radius): break
            if px < 1 or px >= n-1 or py < 1 or py >= n-1: break
            
            if sticking(px, py, grid_ani_2d, stickiness):
                grid_ani_2d[px, py] = 1
                if dist > r_max:
                    r_max = dist
                    spawn = r_max + 10
                break # particle stuck, move to next particle
                
    im2d.set_array(grid_ani_2d)
    ax.set_title(f"2D DLA : Particles={frame}, S={stickiness}")
    return [im2d]

# saves a frame every 50 particles for the 2D triangular DLA animation
def update_tri(frame):
    global spawn, r_max
    particles_per_frame = 50 
    
    for _ in range(particles_per_frame):
        theta = ra.uniform(0, 2*np.pi)
        px, py = int(center + spawn*np.cos(theta)), int(center + spawn*np.sin(theta))
        
        while True:
            px, py = move_2Dt(px, py)
            # kill radius
            dist = np.sqrt((px-center)**2 + (py-center)**2)
            if dist > (spawn * kill_radius): 
                break
                
            if px < 1 or px >= n-1 or py < 1 or py >= n-1: 
                break
            
            if sticking_2Dt(px, py, grid_ani_tri, stickiness, n):
                grid_ani_tri[px, py] = 1
                if dist > r_max:
                    r_max = dist
                    spawn = r_max + 10
                break 

    im_tri.set_array(grid_ani_tri)
    ax.set_title(f"Triangular DLA : Particles={frame}, S={stickiness}")
    return [im_tri]

# saves a frame every 20 particles for the 3D DLA animation
def update_3d(frame):
    global spawn3
    particles_per_frame = 20 
    
    for _ in range(particles_per_frame):
        phi = ra.uniform(0, 2*np.pi)
        cost = ra.uniform(-1, 1)
        theta = np.arccos(cost)
        px = int(center3 + spawn3 * np.sin(theta) * np.cos(phi))
        py = int(center3 + spawn3 * np.sin(theta) * np.sin(phi))
        pz = int(center3 + spawn3 * np.cos(theta))
        
        while True:
            px, py, pz = move_3D(px, py, pz)
            dist = np.sqrt((px-center3)**2 + (py-center3)**2 + (pz-center3)**2)
            if dist > (spawn3 * kill_radius): break
            if px < 1 or px >= n3-1 or py < 1 or py >= n3-1 or pz < 1 or pz >= n3-1: break
            
            if np.any(grid_ani_3d[px-1:px+2, py-1:py+2, pz-1:pz+2] > 0):
                grid_ani_3d[px, py, pz] = 1
                if dist > spawn3 - 2: spawn3 = dist + 5
                break
    
    ax.clear()
    z, y, x = np.nonzero(grid_ani_3d)
    ax.scatter(x, y, z, s=1, c=z, cmap='magma')
    ax.set_title(f"3D DLA : Particles={frame}, S={stickiness}")

# initial variables
n = 250
grid_2d = np.zeros((n, n))
center = n // 2
grid_2d[center, center] = 1
spawn, r_max = 10, 0
num_particles = 5000
stickiness = 0.6
kill_radius = 4

# ####################################################################################
# SECTION 2: DLA 2D PLOT
# Have the particles move, stick, then plot the resulting fractal.
# ####################################################################################

print("")
print("Generating 2D DLA Plot...")

for i in range(num_particles):
    theta = ra.uniform(0, 2*np.pi)
    px, py = int(center + spawn*np.cos(theta)), int(center + spawn*np.sin(theta))
    while True:
        px, py = move_2D(px, py)
        dist = np.sqrt((px - center)**2 + (py - center)**2)
        
        # if too far, stop this particle
        if dist > (spawn * kill_radius):
            break    
        if px < 1 or px >= n-1 or py < 1 or py >= n-1: break
        if sticking(px, py, grid_2d, stickiness):
            grid_2d[px, py] = 1
            dist = np.sqrt((px-center)**2 + (py-center)**2)
            if dist > r_max:
                r_max = dist
                spawn = r_max + 10
            break

print("")
plt.figure(figsize=(6,6))
plt.imshow(grid_2d, cmap='magma')
plt.title(f"2D DLA : Particles={num_particles}, S={stickiness}")
plt.savefig("dla_2d.png", dpi=300)
plt.show()

# ####################################################################################
# SECTION 3: DLA 2D ANIMATION
# Save the growing fractal as frames in an animation.
# ####################################################################################

print("Generating 2D DLA Animation...")
grid_ani_2d = np.zeros((n, n))
grid_ani_2d[center, center] = 1
spawn, r_max = 10, 0

fig, ax = plt.subplots(figsize=(6,6))
im2d = ax.imshow(grid_ani_2d, cmap='magma', animated=True)

save_on_this_particle = 50
ani2d = FuncAnimation(fig, update_2d, frames=range(0, num_particles, save_on_this_particle), interval=50, blit=False)
ani2d.save('dla_2d_animation.gif', writer='pillow')
plt.close()

# ####################################################################################
# SECTION 4: DLA TRIANGULAR 2D PLOT
# ####################################################################################

print("")
print("Generating Triangular DLA Plot...")
grid_tri = np.zeros((n, n))
grid_tri[center, center] = 1
spawn, r_max = 10, 0

for i in range(num_particles):
    theta = ra.uniform(0, 2*np.pi)
    px, py = int(center + spawn*np.cos(theta)), int(center + spawn*np.sin(theta))
    while True:
        px, py = move_2Dt(px, py)
        if px < 1 or px >= n-1 or py < 1 or py >= n-1: break
        if sticking_2Dt(px, py, grid_tri, stickiness, n):
            grid_tri[px, py] = 1
            dist = np.sqrt((px-center)**2 + (py-center)**2)
            if dist > r_max:
                r_max = dist
                spawn = r_max + 10
            break

print("")
plt.figure(figsize=(6,6))
plt.imshow(grid_tri, cmap='magma')
plt.title(f"Triangular DLA : Particles={num_particles}, S={stickiness}")
plt.savefig("dla_triangular.png", dpi=300)
plt.show()

# ####################################################################################
# SECTION 5: DLA TRIANGULAR 2D ANIMATION
# ####################################################################################

print("")
print("Generating Triangular DLA Animation...")
grid_ani_tri = np.zeros((n, n))
grid_ani_tri[center, center] = 1
spawn, r_max = 10, 0

fig, ax = plt.subplots(figsize=(6,6))
im_tri = ax.imshow(grid_ani_tri, cmap='magma', animated=True)

# Update the animation call to use the sub-sampling range
ani_tri = FuncAnimation(fig, update_tri, frames=range(0, num_particles, 50), interval=50, blit=False)
ani_tri.save('dla_triangular_animation.gif', writer='pillow')
plt.close()

# ####################################################################################
# SECTION 6: DLA 3D PLOT
# ####################################################################################

# sections 6 and 7 will switch to 1000 particles to keep the run time reasonable

print("")
print("Generating 3D DLA Plot...")
n3 = 80
grid_3d = np.zeros((n3, n3, n3))
center3 = n3 // 2
grid_3d[center3, center3, center3] = 1
spawn3 = 5

for i in range(num_particles - 4000):
    phi = ra.uniform(0, 2*np.pi)
    cost = ra.uniform(-1, 1)
    theta = np.arccos(cost)
    px = int(center3 + spawn3 * np.sin(theta) * np.cos(phi))
    py = int(center3 + spawn3 * np.sin(theta) * np.sin(phi))
    pz = int(center3 + spawn3 * np.cos(theta))
    
    while True:
        px, py, pz = move_3D(px, py, pz)
        dist = np.sqrt((px-center3)**2 + (py-center3)**2 + (pz-center3)**2)
        if dist > (spawn3 * kill_radius):
            break
        if px < 1 or px >= n3-1 or py < 1 or py >= n3-1 or pz < 1 or pz >= n3-1: break
        if np.any(grid_3d[px-1:px+2, py-1:py+2, pz-1:pz+2] > 0):
            grid_3d[px, py, pz] = 1
            dist = np.sqrt((px-center3)**2 + (py-center3)**2 + (pz-center3)**2)
            if dist > spawn3 - 2: spawn3 = dist + 5
            break

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
z, y, x = np.nonzero(grid_3d)
ax.scatter(x, y, z, s=1, c=z, cmap='magma')
plt.title(f"3D DLA : Particles={num_particles - 4000}, S={stickiness}")
plt.savefig("dla_3d.png", dpi=300)
plt.show()

# ####################################################################################
# SECTION 7: DLA 3D ANIMATION
# ####################################################################################

print("Generating 3D DLA Animation...")
grid_ani_3d = np.zeros((n3, n3, n3))
grid_ani_3d[center3, center3, center3] = 1
spawn3 = 5

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

save_on_this_particle_3D = 20 
ani3d = FuncAnimation(fig, update_3d, frames=range(0, num_particles - 4000, save_on_this_particle_3D), interval=100, blit=False)
ani3d.save('dla_3d_animation.gif', writer='pillow')
plt.close()

# ####################################################################################
# SECTION 8: CAPACITY DIMENSIONS
# ####################################################################################

print("")
print("Final Capacity Dimensions")
print("-" * 30)
dimension_2d = calculate_capacity_dimension(grid_2d)
dimension_tri = calculate_capacity_dimension(grid_tri)
dimension_3d = calculate_capacity_dimension(grid_3d)

print(f"Capacity Dimension (2D):    {dimension_2d:.4f}")
print(f"Capacity Dimension (2D Triangular):  {dimension_tri:.4f}")
print(f"Capacity Dimension (3D):    {dimension_3d:.4f}")
print("-" * 30)
print("")

# ####################################################################################
# END PROJECT 4
# ####################################################################################