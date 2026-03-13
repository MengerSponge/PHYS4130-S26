---
meta:
    author: Alec Snodgrass
    topic:  ODE Project
    course: TN Tech PHYS 4130
    term:   Spring 2026
---
# ODE Solver Project Writeup

## Introduction
Motivate ODE Solvers


##  Solver Algorithms
Explain how each of these solveers work. Explain, in terms of formula's and math, what the algorithm accomplishes. 
Do not need a full derivation, or theoretical explination, but should explain what is going on. 
### Semplectic Integrator: Verlet
### Runge-Kutta 4(5)
### SciPy Integrator: DOP853


## Phase Space Trajectory and Energy Evolution of a Simple Harmonic Oscillator
Insert images of the phase space plots for each of the solvers. These pictures will help demonstrate their capability and a use case. 
### Phase Space Undamped SHM
### Phase Space Damped SHM
### Energy Evolution Undamped SHM
### Energy Evolution Damped SHM


## Strengths and Weaknesses of Each Solver
This section is dedicated to comparing the different solvers. I will lower the number of points each one gets until there is a difference. I will time the functions to find a difference. The two from scipy are very good so it will take a lot to find where their limits are. 

Magic Commands: %timeit

Plots of the phase/energy for different solvers with LOW nts until there is a difference

Take the nts down and plot the same function for the different nts to compare it

There will be images of the plots that - hopefully- show the difference between the solvers. If it is not in the code yet, it will be simple to add. I will only need to write the plotting code and generate a png to include in the writeup.


## Extensions


##  Addressing Questions
### Attribution
What resources did you use on this assignment? People, websites, books, etc.

### Timekeeping
How long did you spend on this assignment? If you didn't keep an accurate log, an estimate is fine.

### Languages, Libraries, Lessons Learned
 2. What libraries did you use in your submission? Were any of them remarkable? Great to use, super annoying to use, etc?

> [!NOTE]
> This section probably shouldn't more than a few sentences long. Record what you learned and move on!