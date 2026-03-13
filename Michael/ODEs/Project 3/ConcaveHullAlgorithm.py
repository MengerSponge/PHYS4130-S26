import numpy as np
import matplotlib.pyplot as plt

#This progam impliments the kth nearest neighbors algoirthm to determine the
#conacve hull for a list of points.

k = 3

def d(P1,P2): #distance function
    return float(np.sqrt((P1[0] - P2[0])**2 + (P1[1] - P2[1])**2))

def P1_to_P2(P1,P2): #vector from P1 to P2
    return[P2[0] - P1[0], P2[1] - P1[1]]

def angle(P1,P2):
    V = P1_to_P2(P1,P2)
    return np.arctan(V[1]/(-V[0]))

def SmallestY(points): #finds the index of point with the smallest y value
    Y = points[0][1]
    I = 0
    for i in range(0,len(points)):
        if Y > points[i][1]:
            Y = points[i][1]
            I = i
    return I

def NearestNeighbors(Points,point_Index): #this function find the k nearest neighbors of a point
    k_nearest_neighbors = list(range(0,k)) #list of indices of nearest neighbors
    distances = [] #stores the index of a point and its distance to the point of interest

    for i in range(0,len(Points)): #compute the distance to all other point
        if i != point_Index: #skip the original point
            distances.append([i,d(Points[i],Points[point_Index])])

    distances.sort(key = lambda x: x[1]) #sort based of distances

    for i in range(0,k):
        k_nearest_neighbors[i] = distances[i][0] #store the indices

    return k_nearest_neighbors #return the indices

def FarthestRightTurn(points,point_index): #find the point the farthest right turn from the desired point

    return

def ConcaveHull(points):
    Hull = [] #this will store the indices of elements in points that make up the concave hull

    Non_Added_Indices = list(range(0,len(points)))
    index0 = SmallestY(points) #start by finding the lowest point

    Hull.append(index0)
    Non_Added_Indices.remove(index0) #removes this index from the list of indices

    return

points = [[1,2],[4,1],[3,5],[3,2],[5,4],[3,1],[2,2]]
index = SmallestY(points)
neighborsIndex = NearestNeighbors(points,index)
neighbors = [points[j] for j in neighborsIndex]

x, y = zip(*points)
plt.plot(x,y, 'o', label = 'All Points')

x,y = points[index]
plt.plot(x,y,'o', label = 'Points')

x, y = zip(*neighbors)
plt.plot(x,y, 'o', label = 'Neighbors')
plt.show()

