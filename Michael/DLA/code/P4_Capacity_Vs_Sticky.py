import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from P4_Functions import AddPoint, box_counting_dimension

# DLA Simulation =======================================================================================
#simul
# simulation parameters
L = 600   #number of slots on one edge of the array
N = 5000  #number of points to add to the aggregate
S_list = [1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.25, 0.2, 0.15, 0.1, 0.05, 0.01]    #sticknynesss

#Computation parameters anf histogram
color_offset = N/4
center = (L-1)//2

thetaMin = -np.pi
thetaMax = np.pi
numBins = 15

bin_edges = np.linspace(thetaMin, thetaMax, 1+numBins)
angle_histo = np.zeros(numBins, dtype=int)

Capacity_Dimensions = []
for S in S_list:
    print( "Stickyness = ", S)

    r_eff = 2 #inital effective radius
    angle_histo = np.zeros(numBins, dtype=int)
    Aggregate = np.zeros((L,L)) 
    Aggregate[center, center] = 10000 #seed point

    for i in range(N):
        r_eff = AddPoint(Aggregate, r_eff, S, L, angle_histo, bin_edges, i+1, color_offset)
        if (i + 1)%500 == 0:
            print(i + 1, " of ", N, " points. Effective Radius: ", r_eff)
    
    D = box_counting_dimension(Aggregate)[0]
    print("Capacity Dimension = " , D)
    print()
    Capacity_Dimensions.append( D)


plt.plot(S_list, Capacity_Dimensions)
plt.title("Capacity Dimension Vs Sticknyness")
plt.xlabel("Stickyness")
plt.ylabel("Capacity Dimension")
plt.savefig("Capacity_Vs_Sticky.png", dpi = 400)
plt.show()


