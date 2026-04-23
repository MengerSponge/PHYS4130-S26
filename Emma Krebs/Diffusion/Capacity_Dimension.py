import Diffusion_Main
import matplotlib.pyplot as plt

data = [(1, 5000), (0.5, 5000), (0.25, 5000), (0.1, 5000), (0.01, 5000)]
D_array = []
size_array = []
N_array = []

for point in data:
    D, size, N = Diffusion_Main.main(*point)
    D_array.append(D)
    size_array.append(size)
    N_array.append(N)

fig2, ax = plt.subplots()

for s, n, i in zip(size_array, N_array, range(len(data))):
    plt.plot(s, n, label=f'{data[i][0]}')
    plt.legend()
    plt.title("Capacity Dimension (Log-Log)")
    plt.xlabel("Sizes of Boxes")
    plt.ylabel("Number of Particles")
    print(D_array[i][0])

plt.show()

fig3, ax = plt.subplots()

stickiness = [t[0] for t in data]
slopes = [d[0] for d in D_array]

plt.plot(stickiness, slopes)
plt.title("Slope vs Stickyiness")
plt.xlabel("Stickiness (Percentage)")
plt.ylabel("Slope (log-log)")

plt.plot()
