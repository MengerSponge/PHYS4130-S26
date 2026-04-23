import Diffusion_Functions
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.animation import PillowWriter
import os

def main(p, num):

    num_particles = num # Note: Only use even numbers or else you will get error because of your division setup
    grid_value = int(num_particles/2)
    grid_size = [grid_value, grid_value]
    current_maximum = 1
    generation_distance = 2
    center = [grid_size[0]/2, grid_size[1]/2]
    probability = p
    kill_distance = 10
    stuck_locations = []
    colors = []

    # Create our quick array. This will be holding the information of particles nearby
    grid_array = np.zeros((grid_size[0], grid_size[1]), dtype=int)
    count = 0

    # Set initial seed 
    grid_array[grid_size[0]//2][grid_size[1]//2] = 1

    while count < num_particles:

        location = Diffusion_Functions.generation_sphere(current_maximum, generation_distance, center)
        particle = Diffusion_Functions.Particle(location, probability)

        while particle.stuck == False:

            particle.random_walk()
            if Diffusion_Functions.particle_from_center(particle, center) > kill_distance:
                break  # kill particle

            neighbors = Diffusion_Functions.get_neighbors(particle.location)
            neighbor_count = 0
            touching = False

            for point in neighbors:
                # Check if anything is touching 
                if grid_array[point[0]][point[1]] == 1:
                    touching = True
                    neighbor_count += 1

            if touching and neighbor_count == 1: # Avoid overfilling and focus on tip growth by saying neighbor count = 1
                result = particle.sticky()

                if result == True:
                    stuck_locations.append(particle.location)
                    colors.append(count)

                    grid_array[particle.location[0]][particle.location[1]] = 1
                    particle.stuck = True

                    if Diffusion_Functions.particle_from_center(particle, center) > current_maximum:
                        current_maximum = Diffusion_Functions.particle_from_center(particle, center)
                        kill_distance = current_maximum + 20
                    count += 1
                    if count % 20 == 0: # Every 20 counts print the count number
                        print(f"For {p}, Count: {count}, Current Radius: {current_maximum}")
                else: 
                    continue

        del particle

    fig, ax = plt.subplots()

    data = np.array(stuck_locations)
    x, y = data.T
    frames =[]

    def update(i):

        ax.clear()

        scale = 20 + i*0.5 # Scale changes with increasing iterations
        size = max(0.05, 100 / scale)
        
        scatter = ax.scatter(x[:i*10], y[:i*10], c=colors[:i*10], marker='s', s=size, alpha=0.5)
        fig.patch.set_facecolor('black')
        ax.set_facecolor('black')
        ax.tick_params(colors='white')
        # ax.axis('off')
        center_x, center_y = center

        ax.set_xlim(center_x - scale, center_x + scale)
        ax.set_ylim(center_y - scale, center_y + scale)
        ax.set_aspect('equal')

        return scatter

    script_dir = os.path.dirname(os.path.abspath(__file__))

    gif_path = os.path.join(script_dir, f"animation_{p}.gif")
    png_path = os.path.join(script_dir, f"last_frame{p}.png")

    frames = len(x) // 10

    ani = animation.FuncAnimation(fig, update, frames=frames, interval=50)

    ani.save(gif_path, writer=PillowWriter(fps=20), savefig_kwargs={"facecolor": "black"})

    update(frames) # Grab last frame
    ax.set_xlim(center[0] - current_maximum, center[0] + current_maximum)
    ax.set_ylim(center[1] - current_maximum, center[1] + current_maximum)
    fig.savefig(png_path)

    D, size, N = Diffusion_Functions.capacity_dimension(stuck_locations)

    return D, size, N, stuck_locations, center
