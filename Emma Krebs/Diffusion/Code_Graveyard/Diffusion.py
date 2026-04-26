import random
import numpy as np

class Node:
    '''
        Class: Node
        Description: Object class Node that will be used to subdivided the 3D space to make generation time much faster.
        Each node contains information on the center of the subdivision, the range of the grid (how big it is), and 
        an empty array of children. After it is subdivided, this children array will be updated with more nodes.
    '''

    def __init__(self, center, grid, depth, parent):
        self.center = center
        self.grid = grid
        self.depth = depth
        self.leaf = True
        self.parent = parent
        self.children = [None, None, None, None, None, None, None, None]
        self.particles = [] # Each Node will only hold one particle at a time for its leaves

    
class Particle:
    '''
        Class: Particle
        Description: Particle object that stores location, age, and probability. Can call a random walk on itself
        and store the data.
    '''
    def __init__(self, location, probability):
        self.location = location
        self.probability = probability
        self.parent = None
        self.stuck = False

    def random_walk(self):
        numbers = [-1, 0, 1] # Random walk options
        random_walked = random.choices(numbers , k=3)

        # Update the location due to random walk
        for i in range(len(random_walked)):
            self.location[i] += random_walked[i]


def insert_particle(node, particle, n):
    '''
        Either insert the particle into its node or into a child of that node.

        Args:
            node (object): The node that is being considered for insertion.
            particle (object): Particle object from the node that needs to be inserted.
            n (int): Max depth value dictated by our grid size.

        Returns:
            Nothing: Simply updates where the particles should be.
    '''
    # If node is a leaf, 
    if node.leaf:
        if node.depth == 1: # initialization of root
            subdivide(node, n)
            index = get_child_index(node, particle)
            insert_particle(node.children[index], particle, n)
            return
        if node.depth == n: # if this is already at max depth, just append particle. This node holds a longer lsit.
            node.particles.append(particle)
            return
        if len(node.particles) == 0:
            node.particles.append(particle)
            return
        
        if len(node.particles) == 1:
            
            old_particles = node.particles
            node.particles = []

            subdivide(node, n)
            value = get_child_index(node, particle)

            # Rerun insert particle so these two lists can find a new home node
            for old in old_particles:
                insert_particle(node.children[value], old, n)
            insert_particle(node.children[value], particle, n)

    
    # Node is an inner node. It has leaf children and we need to sort it first and recall insert.
    else:
        # Calculate the child node this particle is supposed to go into 
        get_child_index(node, particle)
        index = get_child_index(node, particle)
        insert_particle(node.children[index], particle, n)


def subdivide(node, n):
    '''
        Takes a node and subdivides the grid into 8 more parts to create an octree. Quickens computing
        time by limiting the number of possible neighbors a moving particle has to check for sticking. 

        Args:
            node (class): Object for center and grid of a region of our 3D space
            n (int): Cutoff value for a nodes depth depending on starting grid size.

        Returns:
            Nothing: Updates the nodes directly. Doesn't return anything. 
    '''

    if node.depth < n:
        node.leaf = False # It is no longer a leaf case. It has children now!!!
        depth = node.depth + 1 # Update depth
        grid_value = node.grid[0] / 2
        grid = [grid_value, grid_value, grid_value]
        
        offset_center = [-grid_value/2, grid_value/2] # Need to offset center. Can either add or subtract.

        child_index = 0
        for dx in offset_center:
            for dy in offset_center:
                for dz in offset_center:
                    center_x = node.center[0] + dx
                    center_y = node.center[1] + dy
                    center_z = node.center[2] + dz

                    node.children[child_index] = Node([center_x, center_y, center_z], grid, depth, node)
                    child_index += 1
    print("Children created")


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
    v = random.uniform(0,1)

    # From a source in document to create a sphere with equally likely random points on it
    theta = 2*np.pi*u
    phi = np.arccos(1 - 2*v)
    px = int(np.sin(phi) * np.cos(theta) * radius + center[0])
    py = int(np.sin(phi) * np.sin(theta) * radius + center[1])
    pz = int(np.cos(phi) * radius + center[2])

    location = [px, py, pz]

    return location


def get_child_index(node, particle):
    '''
        Find a child node index for a given particle.

        Args:
            node (object): Parent node
            particle (object): Particle we are looking to place
        
        Returns:
            value (int): child index value
    '''
    
    value = 0
    # Calculate the child node this particle is supposed to go into 
    if particle.location[0] >= node.center[0]:
        value |= 4
    if particle.location[1] >= node.center[1]:
        value |= 2
    if particle.location[2] >= node.center[2]:
        value |= 1
    
    return value
    

def find_node(root, point):
    '''
        Given a point this finds out with leaf node it is in and returns that node.

        Args:
            root (Object): Root node that is connected to all other nodes. Can transverse to find leaves.
            point (Array): Location of neighborhood. Want to check if there's any particles around.
        
        Returns:
            node (Object): Returns node of a certain location
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

        if  point[0] >= node.center[0]:
            value |= 4
        if point[1] >= node.center[1]:
            value |= 2
        if point[2] >= node.center[2]:
            value |= 1
        
        # Example: Say we said yes to all three if statements. Then we have 7 and that represents
        # our positive quadrant for this center. 

        found_node = node.children[value]

        if found_node is None:
            return node, value
        node = found_node

    return node, value


def get_neighbors(root, node):
    '''
        Grabs the surrounding neighbors. Includes diagonals as well for a total of 26.
    '''
    neighborhood = []
    # Sharing the same parent
    parent = node.parent
    if parent is not None:
        for child in parent.children:
            if child is not None and child is not node:
                neighborhood.append(child)

    # Next door nodes
    offsets = [-node.grid[0], 0, node.grid[0]]

    for i in offsets: # For x coordinate
        for j in offsets: # For y coordinate
            for k in offsets: # for z coordinate
                if i == j == k == 0: # Skip this one because its where the particle is
                    continue
                else:
                    grid_x = node.center[0] + i
                    grid_y = node.center[1] + j
                    grid_z = node.center[2] + k

                    location = [grid_x, grid_y, grid_z]

                    found_node, _ = find_node(root, location)
                    if found_node is not None and found_node is not node:
                        if found_node not in neighborhood:
                            neighborhood.append(found_node)

    return neighborhood


def stickiness(location_list, particle):
    '''
        The function that determines if a particle will stick. 

        Args:
            location_list (Array): Array of all neighborhood particle locations
            particle (Object): The object of our moving particle.

        Returns:
            Boolean: Returns true if particle sticks. Returns false if it doesn't.
            Array: Array of the location of connection for this particle to stick to.

    '''
    prob = particle.probability
    # Find all the particles that are close enough to stick
    px, py, pz = particle.location
    for lx, ly, lz in location_list:
        dx = lx - px
        dy = ly - py
        dz = lz - pz

        if dx*dx + dy*dy + dz*dz <= 1:
            if random.uniform(0, 1) <= prob:
                print("Sticky")
                return True, (lx, ly, lz)
            
    return False, None
    

def particle_from_center(particle, center):
    '''
        Measures particle from center.

        Arg:
            particle (object): Particle being measured
            center (array): Center of grid

        Returns: 
            distance (int): Distance from center
    '''
    location = particle.location

    distance = np.sqrt((location[0] - center[0])**2 + (location[1] - center[1])**2 + (location[2] - center[2])**2)
    return distance
