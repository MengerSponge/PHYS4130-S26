---
meta:
    author: Emma Krebs
    topic: Diffusion Limited Aggregation
    course: TN Tech PHYS4130
    term: Spring 2026
---

# Diffusion

## Introduction to Diffusion

- Introduce topic. Explain what diffusion is and how it works.
- Explain random walks
- Explain aggregation in programs like this
- Show what a 2D image is supposed to look like
[Include image from wiki/web]

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
