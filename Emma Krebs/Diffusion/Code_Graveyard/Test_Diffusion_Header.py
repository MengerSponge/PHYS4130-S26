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
        self.stuck = False

    def random_walk(self):
        numbers = [-1, 0, 1] # Random walk options
        random_walked = random.choices(numbers , k=3)

        # Update the location due to random walk
        for i in range(len(random_walked)):
            self.location[i] += random_walked[i]


def subdivide(node, depth_value, center, n):
    '''
        Takes a node and subdivides the grid into 8 more parts to create an octree. Quickens computing
        time by limiting the number of possible neighbors a moving particle has to check for sticking. 

        Args:
            node (class): Object for center and grid of a region of our 3D space
            depth_value (int): Cutoff value for a nodes depth depending on starting grid size.

        Returns:
            None: Updates the nodes directly. Doesn't return anything. 
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

                    node.children[child_index] = Node([center_x, center_y, center_z], grid, depth)
                    child_index += 1
        
        old_particles = node.particles
        node.particles = []

        for p in old_particles:
            insert_particle(node, p, n)

        for child in node.children:
            subdivide(child, depth_value, center, n)


def insert_particle(node, particle, n):
    '''
        Either insert the particle into its node or into a child of that node.

        Args:
            node (object): The node that is being considered for insertion.
            particle (object): Particle object from the node that needs to be inserted.

        Returns:
            None: Simply updates where the particles should be.
    '''
    value = 0
    if node.leaf:
        if node.depth == n or len(node.particles) == 0:
            node.particles.append(particle)
            return None
        else:
            old_particle = node.particles[0]
            node.particles = []
            subdivide(node, node.depth, node.center, n)
            insert_particle(node, old_particle, n)
            insert_particle(node, particle, n)
            return
        
    else: # in the case insertion
        if particle.location[0] >= node.center[0]:
            value |= 4
        if particle.location[1] >= node.center[1]:
            value |= 2
        if particle.location[2] >= node.center[2]:
            value |= 1

        insert_particle(node.children[value], particle, n)


def find_node(root, point, n):
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

