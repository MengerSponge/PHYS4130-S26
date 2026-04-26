import Test_Diffusion_Header
import matplotlib.pyplot as py


# Starting parameters
n = 3 # Number of iterations of our octree. Depends on our grid space.
num_particles = 100 # Number of particles used in simulation

grid_value = 8**n
grid = [grid_value, grid_value, grid_value]
center = [grid_value/2, grid_value/2, grid_value/2]
count = 0
probability = 1 # probability of something sticking when encountering a particle
current_maximum = 0 # Will be updated as our program exapnds. 
generation_distance = 3 # The distance to the surface of a sphere for generating a new particle
kill_distance = 10 # If particle wanders past this point, it is killed
stuck_particles = [] # The number of stuck particles that appear. Once in order, we can generate them.

# Generate the very first seed and node, known as root node
root = Test_Diffusion_Header.Node(center, grid, 0)
seed_particle = Test_Diffusion_Header.Particle(center, probability)
Test_Diffusion_Header.insert_particle(root, seed_particle, n)
Test_Diffusion_Header.subdivide(root, root.depth, root.center, n)
node, _ = Test_Diffusion_Header.find_node(root, seed_particle.location, n)

for p in node.particles:
    print(f"Location of particle: {p.location}")


