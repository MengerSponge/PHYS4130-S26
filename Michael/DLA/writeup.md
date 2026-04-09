
## Introduction

## Algorithms and Theory
1) Fractals are loosely charecterized by having arbitrary complexity as they are zoomed in.
    A) We can instead quantify their difference from objects of classical geometry by lookin at different measures of
       dimension such as the capacity dimension and the topological dimension.
         1) Capacity dimension, aka box counting dimension, is determined by covering the plot in boxes and repeatedly
            shrinking them and counting how many contain a point of the fractal.
         2) Topological dimension is the number of independent directions you can move locally within a space.

2) Sample code and explanation for how the DLA algorithm works. This is left non-specific for now since the algorithm is still being tweaked.
    A) The DLA aggregate is stored as a 2D numpy array. There are three primary functions that are used to compute it.
      1) The first function looks at an input point and checks for neighbors in the surrounding array.
      2) The next function does a step of a random walk for a point input into the array.
      3) The last function uses the first two. It starts a random point and walks it until it is next to a neighbor and
         then adds it to the array.
    B) FIrst plot showing a sample DLA plot.

## Capacity Dimension of the DLA
1) Sample code for calculating the capacity dimension.
    A) Plot of the capacity dimension as a function of stickiness

## Animations for N = 10^6 and other interesting plots
