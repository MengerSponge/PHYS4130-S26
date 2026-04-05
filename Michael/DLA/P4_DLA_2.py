
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

np.random.seed()
# Rough Idea for the code:
#
# The Aggregate will be stored in a 2D numpy array of 1s and 0s. 
# 1 = point is here, 0 = no point is here
# A new point will start at a random spot in the aggregagte
# First, we check if any of its 4 neighboring spots are in the aggregate
# If one of them is, then update that spot in the array to 1 and start the next point
# If not, then do a random walk on the point and check again

# DLA Functions ====================================================================================
def Neighbors(Point, Agg): #Takes in the current point and the aggregate and searches for neighbor
    #recall that 2D arrays are accessed as Arr[row (y), column (x)]
    Y, X = Point #point is passed in as the list [Y, X] 

    if Agg[Y, X] > 0: #this is for the case where a point spaws in the array
        return True 
    
    # with clever list slicing and min/max funtions, we can avoid a lot of cases.

    x_low = max([X - 1, 0])       # If X - 1 goes out of the array, then go to 0
    x_high = min([X + 1,L - 1])   # If X + 1 goes out of the array, then go to L - 1

    y_low = max([Y - 1, 0])       # If Y - 1 goes out of the array, then go to 0
    y_high =  min([Y + 1,L - 1])  # If Y + 1 goes out of the array, then go to L - 1

    # this sum will only be greater than 0 if there is a neighbor.
    # it goes over the point (X,Y) twice and the neighbors once (when not at an edge)

    if np.sum(Agg[Y, x_low:x_high+1]) + np.sum(Agg[y_low:y_high+1, X]) > 0:
        return True
    return False

def Walk(Point): #This can be refined later.
    Y, X = Point

    direction = np.random.choice([0, 1])

    if direction == 0: #horizontal
        if X == 0: #at the left
            X += 1
        elif X == L -1: #at the right
            X -= 1
        else:
            X += np.random.choice([-1, 1])
        
        return [int(Y),int(X)]
    
    if direction == 1: #vertical
        if Y == 0: # at the top
            Y += 1
        elif Y == L -1: #at the bottom
            Y -= 1
        else:
            Y += np.random.choice([-1, 1])

        return [int(Y),int(X)]
    
def AddPoint(Aggregate, color_num):
    S = 0.75 #stickyness
    l = len(Aggregate[:, 0])

    X = np.random.choice(range(0, l))
    Y = np.random.choice(range(0, l))
    
    while Aggregate[Y, X] > 0: #seed points outside of the aggregate
        X = np.random.choice(range(0, l))
        Y = np.random.choice(range(0, l))
    X = int(X)
    Y = int(Y)

    while Neighbors([Y, X], Aggregate) == False:
        Y, X = Walk([Y, X])

        if Neighbors([Y, X], Aggregate) == True:
            s = np.random.uniform(0,1) #check this against stickyness S = 1 guranteses it sticks
            if s > S: #smaller than S, walk again and check later
                Y, X = Walk([Y, X])

    Aggregate[Y, X] = color_num
    return

# DLA Simulation =======================================================================================
L = 150 #number of slots on one edge of the array
N = 3000 #number of points to add to the aggregate

#Create the aggregate
Aggregate = np.zeros((L,L)) 
Aggregate[(L-1)//2][(L-1)//2] = 1 #seed point

fig, ax = plt.subplots()
im = ax.imshow(Aggregate, cmap='magma', vmin = 0, vmax = 2*(N+50))

def update(frame):
    AddPoint(Aggregate, 2*(frame + 100))
    im.set_array(Aggregate)

    if (frame + 1)%25 == 0:
        print(frame + 1, " of ", N, " points")

    return [im]

ani = FuncAnimation(fig, update, frames = N, interval = 10, blit=False, repeat = False) #repeat=False or else it adds points forever

#live animation
plt.show()

if 1 == 0: #this is WIP, so it is turned off for now.
    print("Live animation done. Moving to compute the saved animation.")
    #saved animation
    Aggregate = np.zeros((L,L)) 
    Aggregate[(L-1)//2][(L-1)//2] = 1 #seed point

    fig, ax = plt.subplots()
    im = ax.imshow(Aggregate)

    ani.save("DLA.gif", writer = "pillow", fps = 20)

