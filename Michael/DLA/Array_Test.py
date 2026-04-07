import numpy as np

#elements accessed as Arr[row, column]
if 1 == 0: 
    Arr = np.array(
        [[1, 2, 3],
        [4, 5, 6],
        [7, 8, 9]]
    )

    #print(Arr)
    #list slicing is half-open. The end is not included
    print(np.sum(Arr[0:1, 2:3]))
    print(Arr[0, 0:3]) #slices horizontally [1 2 3]
    print(Arr[0:3, 0]) #slices vertically   [1 4 7]
    print(Arr[0:2, 0:2]) #2x2 sub matrix    [[1 2] [4 5]]
    print(np.sum(Arr[0:2, 0:2]))

    np.random.seed()
    print(np.random.randint(0, 2))


#testing the neighbors
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
        # If I want the diagonals to be neighbors too, the I can use
        # np.sum(Agg[y_low:y_high+1, x_low:x_high+1]) > 0
        return True
    return False

L = 5 #number of slots on one edge of the array
N = 4 #number of points to add to the aggregate

#Create the aggregate
Arr = np.zeros((L,L)) 
Arr[(L-1)//2][(L-1)//2] = 1 #seed point
Arr[0, 0] = 1
Arr[2, 1] = 1
Arr[1, 2]
print(Arr)

truth = np.array([[Neighbors([i, j], Arr) for j in range(L)] for i in range(L)])
print(truth)
