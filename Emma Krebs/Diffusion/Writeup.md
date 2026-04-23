---
meta:
    author: Emma Krebs
    topic: Diffusion Limited Aggregation
    course: TN Tech PHYS4130
    term: Spring 2026
---

# Diffusion

## Introduction to Diffusion

Diffusion-limited aggregation (DLA) is the process where particles undergo a random walk, clustering together to form aggregations in the form of fractal branches. We define a random walk as a path that consists of the succession of random steps in some arbitrary mathematical space. For example, the main program this report describes is in a 2D space, so a particle could make a single random step along the positive and negative x- and y-axes. Doing this repeatedly gives us our random walk. An example of how a 3D and 2D aggregate looks like is shown below:

![Image]()
![Image]()
*Fig. 1) Left, right, and middle Riemann sums for y= f(x) on [1, 8] with 5 subintervals.*

PUT TWO IMAGES HERE EMMA

It is important to note that DLA programs are not limited to the circular and spherical shapes above. They can adhere to straight lines, create brownian trees, and have many other shapes:
https://arxiv.org/html/2504.13400#:~:text=Initially%20proposed%20by%20Witten%20and,;%20jungblut%20;%20halsey%20;%20halsey2%20.

PUT TWO MORE IMAGES

DLAs are not solely a computational novelty either, but often appear in nature. For example, mineral deposits, fungi, lightning bolts, snowflakes, and even ants biting off wall paint all follow a form of diffusion-limited aggregation.

https://www.deconbatch.com/2019/10/the-poor-mans-dla-diffusion-limited.html

There are three main factors that affect how an aggregate forms: seed particle location, stickiness probability, and what we consider ‘neighbors’. Stickiness probability determines how likely a particle will become a part of the aggregate when encountering a ‘stuck’ particle. High probabilities mean it will most likely stick to the first particle it encounters, meaning we have thinner, elongated branches. Lower probabilities allow more time for the particle to travel deeper into the structure on its random walk, so for the same number of particles these structures tend to be condensed with little to no branches (See GIFs and Images Section for More!). This also means that low stickiness eats up computation time because of the increased iterations of the particles’ random walks. Neighbor definitions can also change our aggregate branch shapes. For a 2D space, if we only consider the vertical and horizontal neighbors around our particle, we can get more rigid, cardinal growth along our axes compared to an eight neighbor approach. This program utilizes all eight neighbors for the 2D space. Finally, a seed particle is what determines the start of our accumulation for the aggregation. It is placed at a location to become our source point for the following particles to stick to, and it doesn’t necessarily have to be directly in the center of the particles’ generation range. 

https://medium.com/nerd-for-tech/neighborhood-connections-and-connected-components-cedf922dd383 (source)

Now that we know what a DLA is and how different factors can affect them, let us create our own!

## Summary of Code

### Overview
-Explain basic overview of code. Maybe include a drawing of how everything is connected
[Drawing of how everything is connected]

### Initial Idea
-Explain where your initial idea to use octree came from (Game design)
-Explain how octrees work in 3D space [Include image of 3D octree]
-Explain the two main classes that make up most of the code: Particle and Node

### Main Body and Functions (Longest section)
-Explain the process of the steps of main.py 
-As you progress, explain the functions as they are used. Some will be quicker to explain than others, but some major functions to highlight are:
  -subdivide
  -find_node
  -insert_particle
  -stickiness
  -random_walk (inside of particle class)
[Include drawings showcasing what each function does (?)] An important one would be the bitwise children node and how the program determines that

### Gifs/Images
-Talk about using intevrals to generate gifs
-WHere all the images and gifs of different particle number and stickiness will be
  -Planning to include at least 3 different number counts for 5 different stickiness levels

### Convergence and Stability
-Answer to if it terminates
-What happens if particle diffuses away from structure (Increases runtime but I kill off these particles)

### Capacity Dimension and Topological Dimension
- Answer to questions

## Extensions (Doing two)

### Generating long DLA Animation for 1e6 Particles
-Talk about improvements made to octree (Only found out later that octrees are good for C++ and theoretically python, but not in practice).

### How does behavior change in 3D
-Program will be initially in 3D, so need to find another source to compare to

### 

## Languages, Libraries, Lessons Learned

-octree


## Timekeeping

4/9/26: 14 hours

## Sources

-Separate source list on google docs to be implemented
