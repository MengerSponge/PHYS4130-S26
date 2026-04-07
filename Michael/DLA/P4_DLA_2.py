
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

np.random.seed()
# Rough Idea for the DLA algorithm:
#
# The Aggregate will be stored in a 2D numpy array of non-zeros and 0s.  bigger number = newer
# not 0 = aggregate point is here, 0 = no point is here
# A new point will start at a random spot in the array
# First, we check if any of its 4 neighboring spots are in the aggregate
# If one of them is, then update that spot in the array to the current point total and start the next point
# If not, then do a random walk step on the point and check again

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

    #Summation patterns for finding neightbors: (sum > 0 implies next to a neighbor)
    #Adjacent points only: np.sum(Agg[Y, x_low:x_high+1]) + np.sum(Agg[y_low:y_high+1, X]) > 0
    #Diagonal points too: np.sum(Agg[y_low:y_high+1, x_low:x_high+1]) > 0

    if np.sum(Agg[y_low:y_high+1, x_low:x_high+1]) > 0:
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

r_eff = [1] # I need r_eff to be a mutable object, and this is the easiest way to do it for this program
def AddPoint(Aggregate, point_num): #numpy arrays are mutable, so we can alter them in functions easily
    S = 0.9 #stickyness
    l = len(Aggregate[:, 0]) #dimension of the square the aggregate lives in
    
    #old r_eff formula: r_eff = np.sqrt(1+point_num)*2 + 4
    r = r_eff[0]
    theta = np.random.uniform(0, 6.28318530718)

    x = l//2 + (0.5+r)*np.cos(theta)
    y = l//2 + (0.5+r)*np.sin(theta)

    #edge cases to keep poits in the array
    if x < 0:
        x = 0
    if x > l -1:
        x = l-1

    if y < 0:
        y = 0
    if y > l-1:
        y = l-1

    X = int(np.round(x))
    Y = int(np.round(y))
    
    while Neighbors([Y, X], Aggregate) == False:
        Y, X = Walk([Y, X])

        if Neighbors([Y, X], Aggregate) == True:
            s = np.random.uniform(0,1) #check this against stickyness S = 1 guranteses it sticks
            if s > S: #smaller than S, walk again and check later
                Y, X = Walk([Y, X])

    d = np.sqrt((X - (l//2))**2 + (Y - (l//2))**2) #taxi distance between the new point and the origin
    if d > r:
        r_eff[0] = d

    Aggregate[Y, X] = point_num + 250
    return

# DLA Simulation =======================================================================================
L = 300 #number of slots on one edge of the array
N = 5000 #number of points to add to the aggregate

#Create the aggregate
Aggregate = np.zeros((L,L)) 
Aggregate[(L-1)//2][(L-1)//2] = 1 #seed point

fig, ax = plt.subplots()
im = ax.imshow(Aggregate, cmap='magma', vmin = 0, vmax = N)

def update(frame):
    AddPoint(Aggregate, frame)
    im.set_array(Aggregate)

    if (frame + 1)%25 == 0:
        print(frame + 1, " of ", N, " points. Effective Radius: ", r_eff[0])

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

