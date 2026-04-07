import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

if 1 ==0:
    t = np.linspace(0, 10, 100)
    y = np.sin(t)


    fig, axis = plt.subplots()

    axis.set_xlim([min(t), max(t)])
    axis.set_ylim([-2, 2])

    animated_plot, = axis.plot([], [])

    def update(frame):
        animated_plot.set_data(t[:frame], y[:frame])
        return animated_plot,

    anim = FuncAnimation(fig, update, frames = len(t), interval = 25)
    plt.show()




frames = [np.random.rand(20, 20) for _ in range(50)]

fig, ax = plt.subplots()
im = ax.imshow(frames[0])

def update(frame):
    im.set_array(frame)
    return [im]

ani = FuncAnimation(fig, update, frames=frames, interval=100, blit=False)

plt.show()
