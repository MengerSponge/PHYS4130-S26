import numpy as np
import random

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
        dx, dy = random.choice([[-1, 0], [0, 1], [0, -1], [1, 0]])
        self.location[0] += dx
        self.location[1] += dy
        

    def sticky(self):
        # This particle has the potential to get stuck. Check all particles nearby. 

        prob = self.probability
        if random.uniform(0, 1) <= prob:
            return True
            
        return False # It was not stuck
            


def get_neighbors(array):
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
                    grid_x = array[0] + i
                    grid_y = array[1] + j
                    grid_z = array[2] + k

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
    v = random.uniform(0,1)

    # From a source in document to create a sphere with equally likely random points on it
    theta = 2*np.pi*u
    phi = np.arccos(1 - 2*v)
    px = round(np.sin(phi) * np.cos(theta) * radius + center[0])
    py = round(np.sin(phi) * np.sin(theta) * radius + center[1])
    pz = round(np.cos(phi) * radius + center[2])

    location = [px, py, pz]

    return location


def particle_from_center(particle, center):
    '''
        Measures particle from center in squared terms.

        Arg:
            particle (object): Particle being measured
            center (array): Center of grid

        Returns: 
            distance (int): Distance from center
    '''
    location = particle.location

    distance = np.sqrt((location[0] - center[0])**2 + (location[1] - center[1])**2 + (location[2] - center[2])**2)
    return distance

