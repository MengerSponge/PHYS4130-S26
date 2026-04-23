

# Introduction
As particles move through a medium such as air or water, they can bounce off other particles which effects their motion. This phenomena is referred to as Brownian motion. Now, suppose a seed exists that incoming particles can stick to and once a particle sticks to the seed, it becomes a part of the seed. This process is called diffusion-limited aggregation (DLA), and it is observed in many physical systems such as crystal formation.

DLA can be modeled using a simulation in which Brownian motion is simulated as a random walk. In other words, DLA is a process in which particles take random walks in a region where there exists a seed to which the particles stick. As more particles stick to an aggregate it forms into a fractal with a certain property that pertains to its geometry. This property is known as  capacity dimension, and it can have non-integer values as opposed to topological dimension. It is defined as Eq. 1 and can be estimated by overlaying space with boxes of side length $\epsilon$ and counting the number of boxes, N, that contain part of the aggregate. This process is repeated for decreasing $\epsilon$ then $\ln N(\epsilon)$ vs. $\ln (\frac{1}{\epsilon})$ is plotted where the slope yields the capacity dimension.

$$
D_C=\lim_{\epsilon\to 0}\frac{\ln N(\epsilon)}{\ln\left(\frac{1}{\epsilon}\right)}
\qquad\text{(1)}
$$
where $D_{C}$ is the capacity dimension, $\epsilon$ is the size of boxes, and N is the number of boxes. Furthermore, the expected value of the capacity dimension of a DLA fractal is 1.7 [1]. Accordingly, the capacity dimension of a simulated aggregate can be compared with this value to evaluate the accuracy of the simulation. Capacity dimension can additionally be used to evaluate the effects of certain aspects of the model on its accuracy. For instance, if the simulated flux of incoming particles is not uniform, the capacity dimension would likely stray away from its expected value for a uniform flux system.

Another useful fractal dimension is the mass-radius fractal dimension. It is defined by
$$
N(r) = r^{D}\qquad(2)
$$
where $r$ is the distance from the center of the aggregate, $N(r)$ is the number of occupied elements within radius $r$, and $D$ is the fractal dimension. This dimension essentially describes the density of the aggregate. It has a maximum value of 2 in 2D space which would result from a filled in circle.

In this project, we model diffusion limited aggregation for a 2D system with a uniform flux of incoming particles from an infinite radius. We evaluate the capacity dimension of our simulated aggregate at different values of stickiness and radius. Then, we compare our estimated value to its expected value. In addition, we evaluate the effect of letting the spawn radius approach the aggregate.

# Procedure 
Our simulation code is described in a step by step manner below containing the main steps to model DLA. We use a boolean array to model 2D space in which each element's truth value corresponds to whether it is occupied or not. Next, we define a function that spawns particles on an approximately circular parameter centered around the aggregate. The spawn probability is uniformly distributed about the circle. Then, a kill radius is implemented to improve runtime efficiency. Lastly, random.getrandbits() is used to implement random walks, and the 3x3 neighborhood around the particle is checked for the seed after each walk. This process is repeated until the particle sticks to the aggregate.

1.  Initializing a 2D space array
```python
space = np.zeros((length, length), dtype=bool)            
space[length//2, length//2] = True                        
```

A boolean array is used to represent discrete 2D space (length x length blocks): An element is true or 1 if occupied and false or 0 if empty. The seed is initialized in the center.

2. Defining a function to spawn particles
```python
def spawn(spawn_radius):
    theta = np.random.uniform(0, 2 * np.pi) 
    x, y = int(spawn_radius * np.cos(theta)), int(spawn_radius * np.sin(theta))
    return x, y
```
The spawn function is called to spawn a particle at a random angle on an approximate circular parameter around the aggregate. We played around with different conditions for the radius of the parameter. Mainly we explored two conditions based on the radius of the aggregate's longest branch, $r_{max}$: radius = $r_{max}$ + 10 and radius = 2 * $r_{max}$.

3. Implementing a kill radius
```python
def kill(i, x, y):
    global length, space, heat, spawn_radius 

    r_i = int(np.sqrt((x - length//2)**2 + (y - length//2)**2)) + 1          

    if r_i > 3 * spawn_radius + 10:
        return True
    else:
        return False
```
The kill function prevents particles from wandering far away from the circle. This tremendously improves runtime as particles are not allowed to wander far away which causes them to take much longer to find the aggregate. If the particle wanders more than 10 blocks outside the spawn radius, it is killed and respawned. Additionally, a kill radius does not alter the accuracy of the model because a particle has the same probability of being killed regardless of where it spawns and the spawn probability is uniformly distributed about the spawn parameter.

4. Implementing random walks & checking neighbors
```python
    # Loop until the particle sticks to seed
    while (True):                           
        # RANDOM WALK ALGORITHM
        direction = random.getrandbits(2)
        
        new_x, new_y = x, y

        if direction == 0:
            new_y += 1
        elif direction == 1:
            new_y -= 1
        elif direction == 2:
            new_x -= 1
        else:
            new_x += 1

        # ensures space is not already occupied
        if 0 < new_x < length-1 and 0 < new_y < length-1 and not space[new_x, new_y]:
            x, y = new_x, new_y

        # checks if particle is too far from the seed and kills it if it is
        if kill(i, x, y):

            # respawns particle at random point on circle                                                         
            x, y = spawn(spawn_radius)

            # aligns circle to be centered at the seed                                                   
            x += length // 2                                                            
            
            y += length //2                                                            

            continue                                                     

        # Check if particle is next to the seed (3x3 neighborhood)

        # only checks if particle is not on the edge to avoid errors                                                                 
        if 0 < x < length-1 and 0 < y < length-1:                               
            if ((space[x-1, y] or space[x,   y-1] or space[x, y+1]
                or space[x+1, y]) and not space[x, y]):                         
                # not space ensures spot is not already occupied 

                # Chance of sticking to the seed is the stickiness factor
                if np.random.rand() < stickiness:                               
                    space[x, y] = seed
                    # tracks the age of the particle that sticks to the seed
                    heat[x, y] = i + 1                                          
                    resizing_square(i, x, y)
                    resizing_circle(i, x, y)
                    break                                                     
```
The random walk algorithm chooses a random direction for the particle to walk. The random.getrandbits() function is used for optimal runtime. Then, the 4 neighboring elements are checked if they are part of the aggregate, and the current position is required to not already be occupied. The particle takes random walks until it reaches the aggregate. Once the particle reaches the aggregate, it has a chance of sticking. And once the particle sticks, the loop breaks.

# Analysis 

Animation: I don't think I can embed a mp4?

-------

Plots... notice they look like fractals...

<p align="center">
  <img src="50k_1.47CD_1.67FD.png" alt="Figure 1" width="45%">
  <img src="20k_1.46CD_1.61FD.png" alt="Figure 2" width="45%">
</p>

<p align="center">
  <em>Left: Aggregate for 50,000 particles at a sticking probability of 1. Its estimated capacity and mass-radius fractal dimensions are 1.47 and 1.67 respectively. Right: Aggregate for 20,000 particles at a sticking probability of 1. Its estimated capacity and mass-radius fractal dimensions are 1.46 and 1.61 respectively.</em>
</p>
--------------------------------------------------------------------

Analysis procedure for CD and FD...

Fig: Capacity Dimension Plot (Maybe Dual plot this with the one below)
<p align="center">
  <img src="50k_1.47CD.png" alt="Figure 1" width="70%">
</p>

<p align="center">
  <em>Capacity dimension estimation for 50,000 particles at a sticking probability of 1. </em>
</p>
---------------------------------------

Fig: Radius Fractal Dimension
<p align="center">
  <img src="50k_1.67FD.png" alt="Figure 1" width="70%">
</p>

<p align="center">
  <em>The radius-mass fractal dimension of the aggregate generating by 50,000 particles at a sticking probability of 1. </em>
</p>


------------------------------------------
Stickiness Analysis... density, runtime...

Fig: S vs Capacity Dimension (Maybe dual plot this with the one below)
<p align="center">
  <img src="CDvS_fitted.png" alt="Figure 1" width="70%">
</p>

<p align="center">
  <em>Sticking probability vs Capacity Dimension. The capacity dimension tends to get larger for lower sticking probability. </em>
</p>

------------------------------------------
Fig: S vs Radius Fractal Dimension

<p align="center">
  <img src="SvFD.png" alt="Figure 1" width="70%">
</p>

<p align="center">
  <em>Sticking probability vs Radius-mass fractal Dimension. The fractal dimension tends to get larger for lower sticking probability. </em>
</p>

-----------------------------------
Maybe a short runtime analysis demonstration
```python
```
# Conclusions

# References
[1]

# Appendix
### Changelog Summary
The code described above took a long time to optimize the runtime for large N. Below, we share the trials and errors that led up to the code above.
#### Periodic Boundary Conditions
The first big thing we implemented was periodic boundary conditions to prevent the particles from wandering too far away from the seed. We implemented this using the mod operator as shown below:
``` python
        direction = np.random.choice(['up', 'down', 'left', 'right'])
        if direction == 'up':
            y = (y + 1) % length
        elif direction == 'down':
            y = (y - 1) % length
        elif direction == 'left':
            x = (x - 1) % length
        elif direction == 'right':
            x = (x + 1) % length
```


#### Resizing Boundary
Initially, we also started with a 100 block square spawn parameter centered around the seed. After running a simulation with 1000? particles we observed that the particles accumulated on the edges as they spawned in which is not realistic as the spawn parameter just represents the income of the particles at that point--- they don't just spawn there out of thin air. So, we increased N to avoid this; however, increasing N resulted in a much longer run time since the particles took longer to walk further to the seed. So, then we implemented the resizing functions to keep the spawn parameter at a balanced distance from the seed. The distance was chosen to be 3 times the radius of the aggregate point furthest from the origin. 

```python
def resizing():
    ##############################################
    # Adaptive resizing
    global length, space, heat                                          # I should probably pass these later but rn Im just testing stuff
    
    r_max = np.max(np.abs(np.argwhere(space) - length//2)) + 1          # max distance from the center seed
    if length < 3 * r_max:                                              # condition to adjust length
        new_length = 3 * r_max
        if new_length % 2 == 0:                                         # ensures new length is odd so seed can remain centered
            new_length += 1

        new_space = np.zeros((new_length, new_length), dtype=bool)      # adjusts space to new length dimensions
        new_heat = np.zeros((new_length, new_length))                   # adjusts heat to new length dimensions

        shift = (new_length - length) // 2                              # shift such that the space expands outwards from center

        crystal_indices = np.argwhere(space)                            # ^ adjusts space such that it adds space outwards from center
        for x, y in crystal_indices:
            new_space[x + shift, y + shift] = True
            new_heat[x + shift, y + shift] = heat[x, y]
        
        length = new_length                                             # updates length
        space = new_space                                               # updates space
        heat = new_heat                                                 # updates heat

    ##############################################
```

#### Circular Spawn Parameter
Yet, we still weren't satisfied with the runtime, so we changed the shape of the spawn parameter to a circle instead of a square. This change allowed us to spawn particles closer to the crystal. To visualize this, imagine two concentric squares where the particles pass into the exterior square in a uniformly distributed manner and take random walks. They will not pass through the interior square in a uniformly distributed manner. Meanwhile, for a pair of concentric circles, if the particles pass into the exterior circle in a uniformly distributed manner then they will pass through the interior circle uniformly. Naturally, we also wrote a resizing function for the circle which was based of the maximum distance of a crystal point from the center plus 10.

```python
def generate_space(radius):
    num_edge_blocks = int(2 * np.pi * radius)                                                               # num edge blocks is the integer value of circumerence
    theta = np.linspace(0, 2 * np.pi, num_edge_blocks, endpoint=False)                                      # generates num_edge_blocks angles evenly spaced around the circle
    x, y = np.round(radius * np.cos(theta)).astype(int), np.round(radius * np.sin(theta)).astype(int)       # finds the x and y positions of the points on the circle
    return x, y

x, y = generate_space(radius)
circle_set = np.unique(np.column_stack((x, y)), axis=0)         # pairs x and y defines the circle (np.unique removes any duplicates due to rounding)
```


#### Spawn Function
The early versions of the spawn function would define a segmented curve and pick a random point on that curve. We later changed this to pick a random angle at the set radius which is faster.
```python
def spawn(radius):
    theta = np.random.uniform(0, 2 * np.pi)                                                                 # generates a random angle between 0 and 2pi
    x, y = int(radius * np.cos(theta)), int(radius * np.sin(theta))
    return x, y
```

#### Checking r_max
Initially, we were checking the max radius by using np.max to sift through the distance of each point on the crystal for every loop of N. We improved the efficiency of this search by storing the initial r_max and every new point we check if its greater than the current r_max, and if so, we update r_max.

#### Random Walk Algorithm
Initially, we were using np.random.choice to choose up, down, left or right. We later changed this to np.random.randint which is faster.

#### Neighborhood Check
Initially we were slicing the neighborhood with
```python
        if (space[x-1:x+2, y-1:y+2].any()):
```
yet we found that 
```python
        if 0 < x < length - 1 and 0 < y < length - 1:                           # only checks if particle is not on the edge to avoid errors (aggregate does not approach edge so this is fine)
            if (space[x-1, y-1] or space[x-1, y] or space[x-1, y+1] or          # checks 3x3 neighborhood for crystal
                space[x,   y-1] or space[x,   y] or space[x,   y+1] or
                space[x+1, y-1] or space[x+1, y] or space[x+1, y+1]):
```
is more efficient. 

#### New Boundary Conditions
We found that making the particles bounce off the walls was more efficient than using periodic boundary conditions. And they are parity symmetric so it doesn't change much. 

```python
        # RANDOM WALK ALGORITHM
        direction = np.random.randint(4)

        # particles bounce off walls instead of periodic BC but its effectively the same.
        if direction == 0:        
            if y < length - 1:
                y += 1
        elif direction == 1:
            if y > 0:
                y -= 1
        elif direction == 2:
            if x > 0:
                x -= 1
        else:
            if x < length - 1:
                x += 1
```
