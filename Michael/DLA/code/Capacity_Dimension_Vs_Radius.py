import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from P4_Functions import AddPoint, box_counting_dimension

# DLA Simulation =======================================================================================
#simul
# simulation parameters
L = 600   #number of slots on one edge of the array
N = 7500  #number of points to add to the aggregate
S = 1    #sticknynesss

#Create the aggregate
Aggregate = np.zeros((L,L)) 
center = (L-1)//2
Aggregate[center, center] = 10000 #seed point

#Computation parameters anf histogram
r_eff = 2 #inital effective radius
r = 10
r_list = []
color_offset = N/4

thetaMin = -np.pi
thetaMax = np.pi
numBins = 20

bin_edges = np.linspace(thetaMin, thetaMax, 1+numBins)
angle_histo = np.zeros(numBins, dtype=int)

Capacity_Dimensions = []

for i in range(N):
    r_eff = AddPoint(Aggregate, r_eff, S, L, angle_histo, bin_edges, i+1, color_offset)

    if r_eff > r:
        r += 5
        r_list.append(r_eff)
        D = box_counting_dimension(Aggregate)[0]
        Capacity_Dimensions.append(D)

    if (i + 1)%500 == 0:
        print(i + 1, " of ", N, " points. Effective Radius: ", r_eff, " Capacity Dimension: ", D)

plt.plot(r_list, Capacity_Dimensions)
plt.title("Capacity Dimension Vs Radius")
plt.xlabel("Radius")
plt.ylabel("Capacity Dimension")
plt.savefig("Capacity_Dimension_Vs_Radius.png")
plt.show()

plt.imshow(Aggregate)

coords = np.argwhere(Aggregate)
ymin, xmin = coords.min(axis=0)
ymax, xmax = coords.max(axis=0)
lower = min([xmin, ymin])
upper = max([xmax, ymax])

plt.ylim(lower, upper)
plt.xlim(lower, upper)
title = "DLA S = " + str(S) + " N = " + str(N)
plt.savefig(title, dpi = 400)


