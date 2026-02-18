---
meta:
    author: Emma Krebs
    topic: Numeric Integration Project
    course: TN Tech PHYS4130
    term: Spring 2026
---

# Numeric Integration

## Summary of Code

### Numeric Integration Methods -- Trapezoid Rule 

There are several methods when approximating a definite integral over a limit through Riemann Sums. This is when rectangles are formed by dividing the integral [a, b] into N subintervals from [x<sub>0</sub>, x<sub>1</sub>] to [x<sub>N-1</sub>, x<sub>N</sub>], where x<sub>0</sub> is a and x<sub>N</sub> is b. This approximates the area under the curve for a function f(x) when mulitplied by the height (or point on f(x)). For increasing subintervals N our approximation tends to become more accurate (Excluding cases with functions with alternating peaks above and below the x-axis since it can widely change between the number of subintervals depending on the points chosen. However, it will eventually approach the correct approximation). 

Although all these methods will get you to your correct approximation, not all methods are equivalent in terms of efficency. This will effect both our timing, making more complex integration's approximations longer when trying to reach a certain precision or error, and computer resources for each loop run.

There are five main methods that we have discussed so far: leftpoint, rightpoint, midpoint methods, trapezoid rule, and Simpson's rule. Leftpoint, rightpoint, and midpoint methods inherently follow their namesake such that their rectangles align with the left, right, and middle points, respectively, of the top edge of the rectangle.

 The error of the leftpoint and rightpoint methods decrease at a similar rate and are the worst efficency approximations out of the five with the error decreasing at a linear rate. The midpoint method and trapezoid rule are slightly better decreasing at a quadratic proportionality. The best of the five is Simpson's rule, which combines the weighted sums of the midpoint method and trapezoid rule to get an error proportionality with respect to the fourth power. However, for this project we worked on coding a trapezoid method to solve the following integral:

```math
I = \int_0^2 \mathrm{d}x\, \sin^2\left(\sqrt{100x}\right)
```

To do this, let us go a little more in depth with the trapezoid rule. We mentioned before that the proportionality of decreasing error for the leftpoint and rightpoint methods are approximately the same. Thus, the point of evaluation should not make much of a difference assuming our width of a subinterval is small. Since the errors should be roughly the same in magnitude, but opposite in sign, we can have a better method of integration by averaging both methods. This is where we get our trapezoid rule. That is, to approximate the integral
<<<<<<< Updated upstream
 
```math
I = \int_a^b f(x)\,dx
```
 
=======

'''math
I = \int_a^b f(x)\,dx\n
'''

>>>>>>> Stashed changes
we can do
 
```math
I_T = \frac{1}{2} \left( I_L + I_R\right) = \sum_{i=0}^{N-1} \frac{1}{2} \left[ f(a + ih) + f(a + ih + h) \right] 
```
 
In my code, I created the following definition

```python

    def trapezoid_function(func, a, b, n):
    
        h = (b - a) / n # Width of our subinterval
        sum = 0
        # Loop for number of subintervals
        for i in range(n):
            sum += (1 / 2) * h * (func(a + i * h) + func(a + i * h + h))
    
        return sum
```

to represent this equation. The definition takes in the parameters func (given function), a (starting point of integral), b (ending point of integral), and n (the number of subintervals). Then, it loops over a range of n subintervals and continously sums using the trapezoid equation and our given parameters to get the approximate answer of the integral for n subintervals. While the trapezoid equation isn't called in the main body, it is called by another very important definition: the approximator.

```python

    def approximator(f, start, end, true_value, sig_fig):

        n = 1 # Starting number of subintervals
        sums = [] # Empty array to keep track of the sum of the function using n integrals
        parameter = False # Parameter used to know when to end loop once we have approximated close to the answer
    
        while(parameter == False):
    
            current_sum = trapezoid_function(f, start, end, n)
            error = error_calc(current_sum, true_value)
            sums.append((n, current_sum, str(error) + '%'))
            
            if error <= sig_fig:
                parameter = True
            else:
                n *= 2
    
        return sums
```

This function definition is what passes the function, a, and b to the trapezoid_function. It also passes the true value of the integral (calculated through some other method) and the precision to which we want to know the accuracy of the approximation. It tracks the number of subintervals, the apporximated sum of that specific subinterval number, and the error from the true answer. The approximator function will also stop once an accuracy to whatever given significant figure is reached. I changed up my method of calculation for this cutoff from notebook 3. I originally used half the final decimal place of wanted degree of accuracy (ie. 10^-4 accuracy meant 0.00005 was passed). However, this didn't always result in the accuracy we were interested in, so I found a nice slideshow from Illinois.edu that gave me the corrected version we see above. It looks to find when the relative error is less than or equal to 10<sup>-n + 1</sup> where n is the decimal significant figure. An example of this is if relative error is 10<sup>-2</sup> then the approximation of x has at most three significant figures. The approximator calls on error_calc to find this relative error and then compare it to the passed sig_fig value. Since we were interested in accuracy to the 10^-6 degree, I passed 0.00005 as our sig_fig value to get a value within that error. Until it finds that n, it will double the number of subintervals each loop. The result of this is the following table:

            Subintervals     Summation   Error(%)
        0              1  0.9999753124   0.56948%
        1              2  0.7959466253  20.85666%
        2              4  0.6983700870  30.55898%
        3              8  1.0349702802   2.91018%
        4             16  0.9467001204   5.86679%
        5             32  0.9784652387   2.70829%
        6             64  0.9979096693   0.77487%
        7            128  1.0036893156   0.20018%
        8            256  1.0051951162   0.05045%
        9            512  1.0055754278   0.01264%
        10          1024  1.0056707479   0.00316%
        11          2048  1.0056945931   0.00079%
        12          4096  1.0057005553    0.0002%
        13          8192  1.0057020459     5e-05%

### Gaussian Quadrature

As we talked about in the first section, we discovered a better integration method through averaging the leftpoint and rightpoint methods. Now, consider Simpson's rule which fits a parabola using three points from a combination of the midpoint and trapezoid methods. This improved our efficency, but we can improve it even more by adjusting the sample points themselves. By mapping our original [a, b] range to [-1, 1], we can seek optimal x-coordinates at which to evaulate the given function through a weight function. This range is typically the best for Legendre and Chebyshev polynomials, which we will talk about in the next section. In order to complete this mapping, we will use the following u-substitution:

```math
u=\frac{2x-a-b}{b-a}.
```
which also impacts our differential. This u-substitution means we need a du such that

```math
du=\frac{2}{b-a}.
```

Both of these equations were written as function definitions such that

```python
'''
    Definition: u (gaussian_quadrature)
    Parameters: x (Value it needs to convert), a (start), and b (end)
    Description: Maps a, b to [-1, 1] interval. By doing so, you can convert an x in that range.
'''
def u(x, a, b):
    return (2*x - a - b) / (b - a)

'''
    Definition: du (gaussian_quadrature)
    Parameters: a (start), and b (end)
    Description: Replaces du so any integral can use dx as the respective variable.
'''
def du(a, b):
    return 2 / (b - a)
```
To make sure the u function definition worked, I randomly selected a range to map onto [-1, 1]. Thus, as an example of what this mapping looks like we have:

Test for Gaussian Quadrature with a = 1, b = 7, and x = 5\
The starting point 1 becomes:  -1.0\
The ending point 7 becomes:  1.0\
The test variable 5 becomes:  0.3333333333333333

### Legendre Polynomials



![Image](CompProject1Figure.png)

Fig. 1) This is an image depicting the plots of two Legendre polynomials along with their product, P(i)*P(j). The subgraphs start at P(1) and P(1) and increase in i and j along the rows and columns to create a 4x4 group of subplots ranging from P(1) to P(4).

The resulting integration from -1 to 1 for this is given in the following statement:

    P(i)*P(j)             Value  Approx. Value
    0      P1*P1  6.6666666667e-01        0.66667
    1      P1*P2  0.0000000000e+00        0.00000
    2      P1*P3 -8.3266726847e-17       -0.00000
    3      P1*P4  3.7253204139e-16        0.00000
    4      P2*P1  0.0000000000e+00        0.00000
    5      P2*P2  4.0000000000e-01        0.40000
    6      P2*P3  0.0000000000e+00        0.00000
    7      P2*P4 -5.0827397846e-16       -0.00000
    8      P3*P1 -8.3266726847e-17       -0.00000
    9      P3*P2  0.0000000000e+00        0.00000
    10     P3*P3  2.8571428571e-01        0.28571
    11     P3*P4  2.5425414289e-16        0.00000
    12     P4*P1  3.7253204139e-16        0.00000
    13     P4*P2 -5.0827397846e-16       -0.00000
    14     P4*P3  2.5425414289e-16        0.00000
    15     P4*P4  2.2222222222e-01        0.22222

Thus, confirming that any integral P(i)*P(j) where i does not equal j is 0. 

This orthogonality property is what will allow us to define the Legendre's roots and weights as points for Gaussian Quadrature. It acts similar to Simpson's rule with how it weights points differently depending on the optimal method of minimizing error. The orthogonality property is important because Legendre polynomials are orthogonal to all polynomials of a lower degree. Using the roots as points means we can minimize this error because it is proportional to our polynomial roots, meaning any polynomial with a lower degree is 0 because of this orthogonality principle. Kind of similar to the step we took to get to the trapezoid rule by minimizing the leftpoint and rightpoint errors. Therefore, we can calculate exact polynomials for degrees of polynomials 2n-1 or less. From Dr. Reid's NumericReps file, we can write the approximation for an integral as the following summation:

```math
\int_{-1}^{1} \mathrm{d}x\, f(x) \approx \sum_{i=1}^N c_{N,i} f\left(x_{N,i}\right)
```

where the points $`x_{N,i}`$ are the roots of the Nth order legendre polynomial, and the weights are given by the following integral:

```math
c_{i,n}=\frac{1}{P_n^{\prime}(x_{N,i})}\int_{-1}^1\frac{P_n(x)}{x-x_{N,i}} \mathrm{d}x
```

As an additional source, the following diagram depicts the difference in weighting for points between the trapezoid, Simpson's, and Gaussian Quadrature rules:


Fig. 2) The interpolation polynomials for various quadrature rules.

An example of some of these roots and weights are given in the following table for the following Legendre polynomials one through four:

        P(x)                                              Roots                                            Weights
    0     1                                              [0.0]                                              [2.0]
    1     2          [-0.5773502691896257, 0.5773502691896257]                                         [1.0, 1.0]
    2     3     [-0.7745966692414834, 0.0, 0.7745966692414834]  [0.5555555555555558, 0.8888888888888883, 0.555...
    3     4  [-0.8611363115940526, -0.3399810435848563, 0.3...  [0.3478548451374538, 0.6521451548625462, 0.652...


### Extension 2: Challenge Boogaloo

I chose to answer the following question:\
Why are the optimal points for an $N$ order Gaussian quadrature the zeros of $P_N$?

We can start 

## Languages, Libraries, Lessons Learned

The primary language for this assignment was Python where we used the libraries scipy, numpy, and pandas. We have consistenly worked in Python from the beginning of this module to the final project. The pandas library was especially useful in creating tables and organizing information. I also learned how to create subplots. I knew you could make a 2x2 grid of plots, but I didn't know you could make sizes up to 4x4, so that was neat. 

## Timekeeping

As of 2/17/26, I have spent ~28 hours or so on this project. 

## Sources

People Used:

Cricket/Cordell recommending how to make 4x4 subplot table look better. Dr. Reid is also another contributor because I used your raw markdown file's syntax as a base for my equations.

Websites Used:

https://stackoverflow.com/questions/455612/limiting-floats-to-two-decimal-points (For round function)
https://www.desmos.com/calculator/d9rmt4wfoa (Checked trap rule against)
https://www.geeksforgeeks.org/python/printing-lists-as-tabular-data-in-python/ (Table building)
https://cs357.cs.illinois.edu/textbook/assets/slides/03-Errors.pdf (New method for error relation to sigfig)
https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.legendre.html (Legendre polynomials)
https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.set_xlabel.html (More on axes and matplotlib)
https://personal.math.ubc.ca/~israel/m210/lesson19.pdf (Double checked my error for the integration methods)

All the following sources are for extension 2:
https://en.wikipedia.org/wiki/Quadrature_(mathematics)
https://people.math.wisc.edu/~chr/am205/notes/am205_gauss_quad.pdf
https://mathworld.wolfram.com/GaussianQuadrature.html
https://aalexan3.math.ncsu.edu/articles/gauss_quad.pdf
https://phys.libretexts.org/Bookshelves/Astronomy__Cosmology/Celestial_Mechanics_(Tatum)/01%3A_Numerical_Methods/1.16%3A_Gaussian_Quadrature_-_Derivation
https://en.wikipedia.org/wiki/Gaussian_quadrature#Gauss%E2%80%93Legendre_quadrature
https://www.researchgate.net/figure/The-interpolation-polynomials-for-various-quadrature-rules-in-Section-5_fig1_338434277

Books Used:

Python All-In-One for dummies by John C. Shovic, PHD and Alan Simpson. (Mostly emotional support...)
