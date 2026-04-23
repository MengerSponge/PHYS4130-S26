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

<img src="Gifs_and_Images/Of7_p0001_15h.jpg" width="400"> <img src="Gifs_and_Images/DLA_3D.jpg" width="600">

*Fig. 1) Left image is of a basic 2D DLA and the right is a 3D DLA.*

It is important to note that DLA programs are not limited to the circular and spherical shapes above. They can adhere to straight lines, create brownian trees, and have many other shapes:

<img src="Gifs_and_Images/Lungs_FLA.png" width="400"> <img src="Gifs_and_Images/Rec8_3kc2p.jpg" width="600">

*Fig. 2) Miscellanous DLA aggregations of different shapes. The first is just a creative project from a user on X and the second is a Brownian Tree from the WIki*

DLAs are not solely a computational novelty either, but often appear in nature. For example, mineral deposits, fungi, lightning bolts, snowflakes, and even ants biting off wall paint all follow a form of diffusion-limited aggregation.

<img src="Gifs_and_Images/Ants_Eat_Paint_YUM.jpg" width="400"> <img src="Gifs_and_Images/DLA_Cluster.jpeg" width="600">

*Fig. 3) Left image is DLA pattern from ants eating paint and right image is grown from a copper sulfate solution in an electrodeposition cell.*

There are three main factors that affect how an aggregate forms: seed particle location, stickiness probability, and what we consider ‘neighbors’. Stickiness probability determines how likely a particle will become a part of the aggregate when encountering a ‘stuck’ particle. High probabilities mean it will most likely stick to the first particle it encounters, meaning we have thinner, elongated branches. Lower probabilities allow more time for the particle to travel deeper into the structure on its random walk, so for the same number of particles these structures tend to be condensed with little to no branches (See GIFs and Images Section for More!). This also means that low stickiness eats up computation time because of the increased iterations of the particles’ random walks. Neighbor definitions can also change our aggregate branch shapes. For a 2D space, if we only consider the vertical and horizontal neighbors around our particle, we can get more rigid, cardinal growth along our axes compared to an eight neighbor approach. This program utilizes all eight neighbors for the 2D space. Finally, a seed particle is what determines the start of our accumulation for the aggregation. It is placed at a location to become our source point for the following particles to stick to, and it doesn’t necessarily have to be directly in the center of the particles’ generation range. 

Now that we know what a DLA is and how different factors can affect them, let us create our own!

## Summary of Code

### Overview

This DLA program is made up of three files: Capacity_Dimension.py, Diffusion_Main.py, and Diffusion_Functions.py. Capacity_Dimension is what dictates how many particles will run and with what probability of stickiness they will adhere to the aggregate. To calculate the capacity dimension (a measure of the number of particles compared to the aggregate’s cluster size, or how many particles are packed into the size of the cluster), its main purpose is to streamline the generation of aggregations for different parameters without looping through the entire Diffusion_Main.py. This is so we can easily compare stickiness probability to the capacity dimension. Diffusion_Main.py, as seen by its namesake, is the main contributor to creating these aggregations. It loops through the number of particles and calls all the necessary functions from Diffusion_Functions.py to determine which particles should be added to our cluster. It also animates and saves a .gif and .png for the aggregations. Finally, DIffusion_Functions.py is where the inner-workings of the diffusion program are stored. It contains all the objects and functions for this program to work. 

The main idea of this code is that it takes a given number of particles and creates an array the size of half the number of particles by half the number of particles. All stuck particles will be given a designation equal to one such that we can track the positions of the aggregate through the 0’s and 1’s in the array. A seed particle is spawned directly in the center of the grid. Then, the program enters a while loop that will continue until the cluster has the given number of particles attached to the aggregate. New objects labeled particles are generated on the surface of a sphere through a function called the generation sphere, defined by a user set distance from the maximum size of the aggregate, and these particles keep track of their own location and whether one is stuck or unstuck through their object attributes. It can also call an object defined function to randomly walk itself and directly update its location. It continues this until either 1) It wanders too far from the aggregation and gets deleted by a defined kill distance, or 2) Has a neighboring ‘stuck’ particle. If it does find a neighbor, it calls another object defined function to determine if it sticks depending on the stickiness probability. If it succeeds, the stuck attribute is updated and we can add one to the total number of stuck particles for our aggregate. A new kill and generation distance is created based on the furthest particle, and the cycle continues. Once the cluster is complete, the growth is animated and saved with the last frame being turned into a .png. The aggregation is now complete, and the capacity dimension can now be found. 

Next, let us take a closer look at some of the more important functions!

### Main Body and Functions (Longest section)


### Gifs/Images

Capacity_Dimension.py tested five different stickiness probabilities for 5000 particles each. Their animation and final growth is shown below:


We can see that as our probability decreases, our branches get thicker and the particles tend to be condensed. We can measure the capacity dimension over radius such that we can see this relationship:




### Capacity Dimension and Topological Dimension

As we mentioned before, the capacity dimension represents the number of particles inside of our aggregate's maximum size. In other words, it describes the density of fractal patterns for our problem. The topological dimension is what we traditionally think of for dimension, so for a 2D space it would be equal to 2. What does the capacity dimension look like in our 2D DLA then? We must consider our stickiness probabilities again. We saw previously for high probability there were sparse, branching aggregations while small probabilities had compact, more circular structures. We can study the global capacity dimension over these probabilities such that:



## Extensions

### How does behavior change in 3D

The 3D program, while not too logically different, has a few changes to its behavior and generation time. Since particles are not limited to only a 2D space to move, computation time significantly increased because they can randomly walk in three different directions. Additionally, there was some minor changes to the arrays and generation sphere. The generation sphere is now an actual sphere instead of just a circle, and considers an additional randomly chosen angle. Then, of course, the arrays were increased to three dimensions. A stickiness probability of 1 with 3000 particles was generated and is shown below:

### 

## Languages, Libraries, Lessons Learned

The main language was python and I used the libraries numpy, matplotlib, and os. I developed my skills with using object and classes in python and creating images and animations for aggregations. In particular, I learned how to use the os to more effectively store the .gifs and .pngs. I also learned how to do octrees, but they unfortunately did not end up working for this particular program (but they might be useful for my project over the summer!).

## Timekeeping

As of 4/23/26: 36 hours

## Sources

People Used:

Michael helped me with the capacity dimension. 

Websites (Images and Information):

https://tonybaloney.github.io/posts/why-is-python-so-slow.html

https://cemrehancavdar.com/2026/03/10/optimization-ladder/ 

https://eli.thegreenplace.net/2018/slow-and-fast-methods-for-generating-random-integers-in-python/ 

https://en.wikipedia.org/wiki/Octree

https://vpython.org/

https://www.gut-wirtz.de/dla/improvements.html#:~:text=Outside%20the%20release%20radius%20we,the%20cluster%20during%20one%20step. 

https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.cKDTree.html 

https://pypi.org/project/pyoctree/ 

https://github.com/jcranch/octrees/blob/master/octrees/octrees.py

https://medium.com/data-science/neighborhood-analysis-kd-trees-and-octrees-for-meshes-and-point-clouds-in-python-19fa96527b77 

https://en.wikipedia.org/wiki/Color_quantization 

https://eisenwave.github.io/voxel-compression-docs/svo/svo.html#:~:text=Best%20Case%20for%20Regular%20Octrees,best%20case%20is%20rarely%20encountered. 

https://delimitry.blogspot.com/2016/02/octree-color-quantizer-in-python.html#:~:text=As%20each%20leaf%20has%20the,Delimitry%20at%204:14%20PM 

https://www.geeksforgeeks.org/dsa/octree-insertion-and-searching/ 

https://www.eskimo.com/~scs/cclass/int/sx4ab.html#:~:text=The%20&%20operator%20performs%20a%20bitwise,exclusive%2DOR%20on%20two%20integers. 

https://vispy.org/api/vispy.scene.visuals.html

https://towardsdatascience.com/neighborhood-analysis-kd-trees-and-octrees-for-meshes-and-point-clouds-in-python-19fa96527b77/ 

https://markjstock.org/dla3d/ 

https://discussions.unity.com/t/octree-subdivision-problem-solved/405500 
https://alea.impa.br/articles/v14/14-15.pdf

https://docs.actian.com/ingres/11.0/index.html#page/QUELRef/Numeric_Data_Types.htm 

https://softologyblog.wordpress.com/category/diffusion-limited-aggregation/#:~:text=Another%20thing%20you%20want%20to,hits%20per%20update%20each%20frame.

https://arxiv.org/html/2504.13400#:~:text=Initially%20proposed%20by%20Witten%20and,;%20jungblut%20;%20halsey%20;%20halsey2%20. 

https://www.deconbatch.com/2019/10/the-poor-mans-dla-diffusion-limited.html 

https://medium.com/nerd-for-tech/neighborhood-connections-and-connected-components-cedf922dd383 

https://mathworld.wolfram.com/CapacityDimension.html 

https://archive.lib.msu.edu/crcmath/math/math/c/c040.htm 

https://www.w3schools.com/python/trypython.asp?filename=demo_matplotlib_pyplot 

https://matplotlib.org/stable/api/_as_gen/matplotlib.animation.FuncAnimation.html 
