import Diffusion_Main
import Diffusion_Functions
import matplotlib.pyplot as plt

data = [(1, 5000), (0.5, 5000), (0.25, 5000), (0.1, 5000), (0.01, 5000)]
D_array = []
size_array = []
N_array = []

fig2, ax2 = plt.subplots()
fig3, ax3 = plt.subplots()
fig4, ax4 = plt.subplots()

for point in data:
    D, size, N, stuck_locations, center = Diffusion_Main.main(*point)
    D_array.append(D)
    
    ax2.plot(size, N, label=f'{point[0]}')
    print(D[0])

    radii, Ds = Diffusion_Functions.capacity_dimension_vs_radius(stuck_locations, center)
    ax3.plot(radii, Ds, marker='o', label=f'{point[0]}')

ax2.set_title("Box-Counting Scaling (Log–Log)")
ax2.set_xlabel("log(1 / box size)")
ax2.set_ylabel("log(N)")
ax2.legend()

ax3.set_xscale('log')
ax3.set_title("Fractal Dimension vs Radius")
ax3.set_xlabel("Radius (r)")
ax3.set_ylabel("D(r)")
ax3.legend()

stickiness = [t[0] for t in data]
slopes = [d[0] for d in D_array]

ax4.plot(stickiness, slopes)
ax4.set_title("Fractal Dimension vs Stickiness")
ax4.set_xlabel("Stickiness (Percentage)")
ax4.set_ylabel("Fractal Dimension")

plt.show()
