
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
from P4_Functions import AddPoint, box_counting_dimension

# DLA Simulation =======================================================================================
# simulation parameters
L = 600   #number of slots on one edge of the array
N = 5000  #number of points to add to the aggregate
S = 1 #sticknynesss

#Create the aggregate
Aggregate = np.zeros((L,L)) 
center = (L-1)//2
Aggregate[center, center] = 10000 #seed point

fig, ax = plt.subplots()
ax.axis("off")
im = ax.imshow(Aggregate, cmap='magma', vmin = 0, vmax = N)

#Computation parameters anf histogram
r_eff = 2 #inital effective radius
color_offset = N/4

thetaMin = -np.pi
thetaMax = np.pi
numBins = 15

bin_edges = np.linspace(thetaMin, thetaMax, 1+numBins)
angle_histo = np.zeros(numBins, dtype=int)

aggList = [np.copy(Aggregate)] #store computed frames

# Live Simulation
def update(frame):
    global r_eff
    r_eff = AddPoint(Aggregate, r_eff, S, L, angle_histo, bin_edges, frame+1, color_offset)
    im.set_array(Aggregate)

    coords = np.argwhere(Aggregate)

    ymin, xmin = coords.min(axis=0)
    ymax, xmax = coords.max(axis=0)
    lower = min([xmin, ymin])
    upper = max([xmax, ymax])

    ax.set_xlim(lower, upper)
    ax.set_ylim(lower, upper)

    if (frame + 1)%25 == 0:
        print(frame + 1, " of ", N, " points. Effective Radius: ", r_eff)
        aggList.append(np.copy(Aggregate))
    return [im]

ani = FuncAnimation(fig, update, frames = N, interval = 10, blit=False, repeat = False) #repeat=False or else it adds points forever

plt.show()

#saving a video of the simulation
fig2, ax2 = plt.subplots()
ax2.axis("off")
im2 = ax2.imshow(aggList[0], cmap='magma', vmin = 0, vmax = N)

def update2(frame):
    im2.set_array(aggList[frame])

    coords = np.argwhere(aggList[frame])

    ymin, xmin = coords.min(axis=0)
    ymax, xmax = coords.max(axis=0)
    lower = min([xmin, ymin])
    upper = max([xmax, ymax])

    ax2.set_xlim(lower, upper)
    ax2.set_ylim(lower, upper)
    return [im2]

ani = FuncAnimation(fig2, update2, frames = int(N/25), interval = 10, blit=False, repeat = False) #repeat=False or else it adds points forever

writer = FFMpegWriter(fps=30, bitrate=1800)
ani.save("animation.mp4", writer=writer)



#last stuff to output
print("Relative frequencies for the angular distribution (should be uniform)")
print(angle_histo)
plt.plot(bin_edges[0:len(bin_edges)-1], angle_histo)
D, sizes, counts = box_counting_dimension(Aggregate)

print("Capacity Dimension")
print(D)