import Diffusion
import matplotlib.pyplot as plt
import numpy as np

# ------- Starting Parameters -------
n = 3 # Number of iterations of our octree. Depends on our grid space.
num_particles = 1000 # Number of particles used in simulation

grid_value = 8**n
grid = [grid_value, grid_value, grid_value]
center = [grid_value/2, grid_value/2, grid_value/2]
probability = 1 # probability of something sticking when encountering a particle
stuck_particles = [] # The number of stuck particles that appear. Once in order, we can generate them.
current_maximum = 1 # Will be updated as our program exapnds. 
generation_distance = 5 # The distance to the surface of a sphere for generating a new particle
kill_distance = 10 # If particle wanders past this point, it is killed
stuck_particles = set() # The number of stuck particles that appear. Once in order, we can generate them.

# Generate the very first seed and root

root = Diffusion.Node([int(center[0]), int(center[1]), int(center[2])], grid, 1, None)
seed_particle = Diffusion.Particle([int(center[0]), int(center[1]), int(center[2])], probability)
Diffusion.insert_particle(root, seed_particle, n)

count = 0
while count < num_particles:

    location = Diffusion.generation_sphere(current_maximum, generation_distance, center)
    particle = Diffusion.Particle(location, probability)

    # Diffusion.insert_particle(root, particle, n)

    while particle.stuck is False:

        particle.random_walk()
        if Diffusion.particle_from_center(particle, center) > kill_distance:
            break  # kill particle

        node, _ = Diffusion.find_node(root, particle.location)

        total_location_set = set() # an unordered collection of unique elements. good for here

        if node is not None:
            for p in node.particles:
                total_location_set.add(tuple(p.location)) # Can only add tuples, not arrays
            neighbors = Diffusion.get_neighbors(root, node)

            # Now visit our neighbors and see if we can get surrounding node's particles
            for neighbor in neighbors:
                for p in neighbor.particles:
                    total_location_set.add(tuple(p.location))
            total_location_list = list(total_location_set)

            test, parent_location = Diffusion.stickiness(total_location_list, particle)
            if test:
                stuck_particles.add(particle)
                particle.parent = parent_location
                particle.stuck = True
                count += 1

                # New particle means we have to check if we can subdivide.
                Diffusion.insert_particle(root, particle, n)
        
data = np.array(stuck_particles)
x, y = data.T
plt.scatter(x, y, marker='s', s=1, alpha=0.5, edgecolors=None)
plt.gca().set_facecolor('black')
plt.axis('equal')
plt.xlim(0, grid[0])                                                                                                       
plt.ylim(0, grid[1])
plt.show()
    



