
import numpy as np
import random
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

np.random.seed()

# DLA Functions ====================================================================================
def Neighbors(Agg, coords, Length): #Takes in the current point and the aggregate and searches for neighbor
    X, Y = coords #coordinates passed in as the list [X, Y] 

    x_low = max([X - 1, 0])       # If X - 1 goes out of the array, then go to 0
    x_high = min([X + 1,Length - 1])   # If X + 1 goes out of the array, then go to L - 1

    y_low = max([Y - 1, 0])       # If Y - 1 goes out of the array, then go to 0
    y_high =  min([Y + 1,Length - 1])  # If Y + 1 goes out of the array, then go to L - 1

    #diagonal neigbors
    #np.count_nonzero(Agg[y_low:y_high+1, x_low:x_high+1])

    #4 neighbors 
    #np.count_nonzero(Agg[y_low:y_high+1, X]) + np.count_nonzero(Agg[Y, x_low:x_high+1]) -  np.count_nonzero(Agg[Y,X])

    #3 neighbors
    #np.count_nonzero(Agg[Y:y_high+1, X]) + np.count_nonzero(Agg[Y, x_low:x_high+1]) -  np.count_nonzero(Agg[Y,X])

    return np.count_nonzero(Agg[y_low:y_high+1, x_low:x_high+1])

def Walk(Agg, coords, Length): #This can be refined later.
    X, Y = coords

    direction = random.choice([[1,1], [1,-1], [-1,1], [-1,-1], [0, 1], [0,-1], [1,0], [-1,0]])
    while Agg[Y+direction[1], X + direction[0]] > 0:
        direction = random.choice([[1,1], [1,-1], [-1,1], [-1,-1], [0, 1], [0,-1], [1,0], [-1,0]])

    X = X + direction[0]
    Y = Y + direction[1]

    if X < 0: #at the left
        X = 0
    if X > Length-1: #at the right
        X =  Length-1

    if Y < 0: # at the top
        Y = 0
    if Y > Length-1: #at the bottom
        Y = Length -1
        
    return [int(X),int(Y)]

def CalcAngle(coords, Length):
    X, Y = coords
    X = X - (Length-1)//2
    Y = Y - (Length-1)//2
    return np.atan2(Y, X)

def ChooseAngle(angle_histo, bin_edges):
    min_count = angle_histo.min()
    candidates = np.where(angle_histo == min_count)[0]

    i = np.random.choice(candidates)
    left = bin_edges[i]
    right = bin_edges[i+1]

    return np.random.uniform(left, right)

def StartPoint(start_radius, Length, angle_histo, bin_edges):

    theta = ChooseAngle(angle_histo, bin_edges)

    x = (Length-1)//2 + (start_radius)*np.cos(theta)
    y = (Length-1)//2 + (start_radius)*np.sin(theta)

    #edge cases to keep poits in the array
    if x < 0:
        x = 0
    if x > Length-1:
        x = Length-1

    if y < 0:
        y = 0
    if y > Length-1:
        y = Length-1

    return [int(x), int(y)]

def AddPoint(Agg, r_eff, Sticky , Length, angle_histo, bin_edges, point_num, color_offset): #numpy arrays are mutable, so we can alter them in functions easily
    r = r_eff

    start_radius = r + 3
    kill_radius = start_radius + 1

    X, Y = StartPoint(start_radius, Length, angle_histo, bin_edges)

    NumNeighbors = Neighbors(Agg, [X,Y], Length)
    while NumNeighbors == 0:
        
        X, Y =  Walk(Agg, [X, Y], Length)

        #checks for new neighbors
        NumNeighbors = Neighbors(Agg, [X,Y], Length)
        if NumNeighbors > 0:
            p_stick = 1 - (1-Sticky)**NumNeighbors
            if p_stick < np.random.uniform(0, 1): #checks to see if the particle stuck
                NumNeighbors = 0

        d = np.sqrt((X - ((Length-1)//2))**2 + (Y - ((Length-1)//2))**2)
        if  d > kill_radius:
            X, Y =  StartPoint(start_radius, Length, angle_histo, bin_edges)

    if d > r:
        r = d
    # add point
    Agg[Y, X] = point_num + color_offset

    #update histogram
    theta = CalcAngle([X, Y], Length)
    bin = 0
    while theta > bin_edges[1+bin]:
        bin += 1
    angle_histo[bin] = angle_histo[bin] + 1

    #return new radius
    return r

def box_counting_dimension(agg, min_box_size=1):

    # Convert to binary mask
    mask = (agg != 0)
    coords = np.argwhere(mask)
    ymin, xmin = coords.min(axis=0)
    ymax, xmax = coords.max(axis=0)

    mask = mask[ymin:ymax+1, xmin:xmax+1]

    N = min(mask.shape)

    # Use powers of 2 for box sizes
    max_power = int(np.log2(N))
    sizes = [2**k for k in range(int(np.log2(min_box_size)), max_power + 1)]

    counts = []

    for size in sizes:
        # Trim array so it divides evenly into size X size boxes
        trimmed = mask[:N - (N % size), :N - (N % size)]

        # Reshape into blocks. It's a 4 element shape that that essentially
        # has the number of blocks along an edge in the 1st and 3rd spot
        # and the length of each box in the 2nd and fourth spot. 
        new_shape = (trimmed.shape[0] // size, size,
                     trimmed.shape[1] // size, size)

        blocks = trimmed.reshape(new_shape)

        # Check if any pixel in each block is occupied
        occupied = blocks.max(axis=(1, 3))

        # Count occupied boxes
        count = np.sum(occupied)
        counts.append(count)

    sizes = np.array(sizes)
    counts = np.array(counts)

    # Remove zeros (log undefined)
    valid = counts > 0
    sizes = sizes[valid]
    counts = counts[valid]

    # Linear fit in log-log space
    coeffs = np.polyfit(np.log(1/sizes), np.log(counts), 1)
    D = coeffs[0]

    return D, sizes, counts
