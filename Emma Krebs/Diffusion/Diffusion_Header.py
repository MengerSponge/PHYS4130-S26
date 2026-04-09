import numpy as np
import random


class Node:
    '''
        Class: Node
        Description: Object class Node that will be used to subdivided the 3D space to make generation time much faster.
        Each node contains information on the center of the subdivision, the range of the grid (how big it is), and 
        an empty array of children. After it is subdivided, this children array will be updated with more nodes.
    '''

    def __init__(self, center, grid, depth):
        self.center = center
        self.grid = grid
        self.depth = depth
        self.leaf = True
        self.children = [None, None, None, None, None, None, None, None]
        self.particles = None # Each Node will only hold one particle at a time for its leaves

    
class Particle:
    '''
        Class: Particle
        Description: Particle object that stores location, age, and probability. Can call a random walk on itself
        and store the data.
    '''
    def __init__(self, location, probability):
        self.location = location
        self.probability = probability
        self.stuck = False

    def random_walk(self):
        numbers = [-1, 0, 1] # Random walk options
        random_walked = random.choices(numbers , k=3)

        # Update the location due to random walk
        for i in range(len(random_walked)):
            self.location[i] += random_walked[i]


def subdivide(node, depth_value, center):
    '''
        Takes a node and subdivides the grid into 8 more parts to create an octree. Quickens computing
        time by limiting the number of possible neighbors a moving particle has to check for sticking. 

        Args:
            node (class): Object for center and grid of a region of our 3D space
            depth_value (int): Cutoff value for a nodes depth depending on starting grid size.

        Returns:
            None: Updates the nodes directly. Doesn't return anything. 
    '''

    if node.depth < depth_value:
        node.leaf = False # It is no longer a leaf case. It has children now!!!
        depth = node.depth + 1 # Update depth
        grid_value = node.grid[0] / 8
        grid = [grid_value, grid_value, grid_value]
        for i in range(8):
            node.children[i] = Node(center/4, grid, depth)
        index = sort_particles_subdivision(node)
        node.children[index].particles = node.particles
        return None
    else:
        return None


def generate_nearby_point(point):
    '''
        Takes in a given point from the moving particle and returns the 26 directions directly surrounding 
        particle so it can search the nodes of these locations. 

        Args:
            point (array of int): Array of coordinates of moving particle 

        Returns:
            point_grid (array of arrays): Returns an array for coordinate points surrounding particle. 

    '''
    point_grid = []

    for i in [-1, 0, 1]: # For x coordinate
        for j in [-1, 0, 1]: # For y coordinate
            for k in [-1, 0, 1]: # for z coordinate
                if i == j == k == 0: # Skip this one because its where the particle is
                    continue
                else:
                    grid_x = point[0] + i
                    grid_y = point[1] + j
                    grid_z = point[2] + k

                    point_grid.append([grid_x, grid_y, grid_z])

    return point_grid


def generation_sphere(current_maximum, generation_distance, center):
    '''
        Generates a location for a new particle based on the current size of DLA and set generation distance
        by the user.

        Args:
            current_maximum (int): The current furthest distance of a particle from the sphere.
            generation_distance (int): Preset value by the user between current_maximum and shell generation
            center (array): Array of center coordintaes

        Returns:
            location (array): Location of where new particle has been spawned in.
    '''

    radius = current_maximum + generation_distance
    u = random.uniform(0,1)

    # From a source in document to create a sphere with equally likely random points on it
    theta = 2*np.pi*u
    phi = np.arccos(1 - 2*u)
    px = round(np.sin(phi) * np.cos(theta) * radius) + center[0]
    py = round(np.sin(phi) * np.sin(theta) * radius) + center[1]
    pz = round(np.cos(phi) * radius) + center[2]

    location = [px, py, pz]

    return location


def kill_or_be_killed(location, kill_radius, max_distance, center):
    '''
        In the case the particle starts to get too far away and might start running up the runtime, it is
        deleted from the program.

        Args:
            location (Array): The three points of the particles location
            kill_radius (int): Defined distance from the max distance that will get a particle killed
            max_distance (int): Furthest particle distance from the center
            center (Array): Center of the original seed particle

        Returns:
            Test (boolean): Tells us if this has exited the region of consideration. 
    '''
    radius_center = np.sqrt((location[0] - center[0])**2 + (location[1] - center[1])**2 + (location[2] - center[2])**2)

    rad = kill_radius + max_distance
    if radius_center >= rad:
        return True
    else:
        return False
    

def find_node(root, point):
    '''
        Given a point this finds out with leaf node it is in and returns that node.

        Args:
            root (Object): Root node that is connected to all other nodes. Can transverse to find leaves.
            point (Array): Location of neighborhood. Want to check if there's any particles around.
        
        Returns:
            locations_of_particle (Array): Returns an array of all nearby particle locations
            value (int): Index of child node for point location
    '''

    # From the wiki, we can use the color quantization program which determines the child node 
    # via the formula 4r + 2g + b, but here instead of red, green, and blue we can use 
    # our postive and negative 3 directions. Thus
    node = root
    value = 0

    # This will loop until it finds the leaf node to extract the particles
    while node.leaf != True:

        value = 0

        if  point[0] >= root.center[0]:
            value |= 4
        if point[1] >= root.center[1]:
            value |= 2
        if point[2] >= root.center[2]:
            value |= 1
        
        # Example: Say we said yes to all three if statements. Then we have 7 and that represents
        # our positive quadrant for this center. 

        found_node = node.children[value]
        if found_node.leaf == True:
            node = found_node
            break
        else:
            node = found_node

    particle_array = node.particles

    location_of_particles = []
    for particle in particle_array:
        location_of_particles.append(particle.location)

    return location_of_particles, value


def stickiness(location_list, particle):
    '''
        The function that determines if a particle will stick. 

        Args:
            location_list (Array): Array of all neighborhood particle locations
            particle (Object): The object of our moving particle.

        Returns:
            Boolean: Returns true if particle sticks. Returns false if it doesn't.

    '''
    # Find all the particles that are close enough to stick
    close_distance = []
    for location in location_list:
        distance = np.sqrt((location[0] - particle.location[0])**2 + (location[1] - particle.location[1])**2 + (location[2] - particle.location[2])**2)
        if distance <= 1:
            close_distance.append(location)

    # If there was nothing close enough, return false
    if len(close_distance) == 0:
        return False 
    else: 
        for close_loc in close_distance:
            prob = particle.probability

            if random.uniform(0, 1) <= prob:
                return True
            else:
                return False


def sort_particles_subdivision(node):
    '''
        Sort particle so its in its new leaf node.

        Args:
            node (object): The node that is getting subdivided and needs to relocate its children
        
        Returns:
            value (int): Index of child node where particle will now reside. 
    '''
    # Need to find index, so we will use color method again.
    particle = node.particles

    location = particle.location

    value = 0

    if location[0] >= node.center[0]:
            value |= 4
    if location[1] >= node.center[1]:
            value |= 2
    if location[2] >= node.center[2]:
            value |= 1

    return value
