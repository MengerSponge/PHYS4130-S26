import Diffusion_Header

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

print("Done here 0")

# Generate the very first seed and node, known as root node
root = Diffusion_Header.Node(center, grid, 0)
seed_particle = Diffusion_Header.Particle(center, probability)
Diffusion_Header.insert_particle(root, seed_particle, n)
Diffusion_Header.subdivide(root, root.depth, root.center, n)
print(Diffusion_Header.find_node(root, seed_particle.location, n))


stuck_particles.append(seed_particle.location)

for i in range(2):

    location = Diffusion_Header.generation_sphere(current_maximum, generation_distance, center)
    particle = Diffusion_Header.Particle(location, probability)
    print(location)

    while particle.stuck != True:
        particle.random_walk()
        nearby_points = Diffusion_Header.generate_nearby_point(particle.location)

        total_location_list = [] # Total list of particles we need to compare with for our point
        for point in nearby_points:
            found_node, _ = Diffusion_Header.find_node(root, point, n)
            print(found_node, found_node.particles, _)

            for p in found_node.particles:
                print("Particle ID:", id(p), "Location:", p.location)
                if p.location not in total_location_list:
                    total_location_list.append(p.location)

        print(total_location_list)

        test = Diffusion_Header.stickiness(total_location_list, particle)
        if test == True: # Particle has stuck
            stuck_particles.append(particle.location)
            particle.stuck = True
            node, _ = Diffusion_Header.find_node(root, particle.location, n)

            # New particle means we have to check if we can subdivide.
            found_node, value = Diffusion_Header.find_node(root, particle.location, n)
            Diffusion_Header.insert_particle(found_node, particle, n)
                        
    if Diffusion_Header.particle_from_center(particle, center) > current_maximum:
        current_maximum = Diffusion_Header.particle_from_center(particle, center)

print(stuck_particles)
