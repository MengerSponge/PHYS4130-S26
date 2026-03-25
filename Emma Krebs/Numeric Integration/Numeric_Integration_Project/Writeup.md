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

There are several methods when approximating a definite integral over a limit through Riemann Sums. This is when rectangles are formed by dividing the integral [a, b] into N subintervals from [x<sub>0</sub>, x<sub>1</sub>] to [x<sub>N-1</sub>, x<sub>N</sub>], where x<sub>0</sub> is a and x<sub>N</sub> is b. This approximates the area under the curve for a function f(x) when multiplied by the height (or point on f(x)). For increasing subintervals N our approximation tends to become more accurate (Excluding cases with functions with alternating peaks above and below the x-axis since it can widely change between the number of subintervals depending on the points chosen. However, it will eventually approach the correct approximation). 

Although all these methods will get you to your correct approximation, not all methods are equivalent in terms of efficiency. This will affect both our timing, making more complex integration's approximations longer when trying to reach a certain precision or error, and computer resources for each loop run.

There are five main methods that we have discussed so far: leftpoint, rightpoint, midpoint methods, trapezoid rule, and Simpson's rule. Leftpoint, rightpoint, and midpoint methods inherently follow their namesake such that their rectangles align with the left, right, and middle points, respectively, of the top edge of the rectangle. An example of the three are show in the following figure.

![Image](LeftRightMid.png)
*Fig. 1) Left, right, and middle Riemann sums for y= f(x) on [1, 8] with 5 subintervals.*

The error of the leftpoint and rightpoint methods decrease at a similar rate and are the worst efficiency approximations out of the five with the error decreasing at a linear rate. The midpoint method and trapezoid rule are slightly better, decreasing at a quadratic proportionality. The best of the five is Simpson's rule, which combines the weighted sums of the midpoint method and trapezoid rule to get an error proportionality with respect to the fourth power. However, for this project we worked on coding a trapezoid method to solve the following integral:

```math
I = \int_0^2 \mathrm{d}x\, \sin^2\left(\sqrt{100x}\right)
```

To do this, let us go a little more in depth with the trapezoid rule. We mentioned before that the proportionality of decreasing error for the leftpoint and rightpoint methods are approximately the same. Thus, the point of evaluation should not make much of a difference assuming our width of a subinterval is small. Since the errors should be roughly the same in magnitude, but opposite in sign, we can have a better method of integration by averaging both methods. This is where we get our trapezoid rule. That is, to approximate the integral
 
```math
I = \int_a^b f(x)\,dx
```
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

to represent this equation. The definition takes in the parameters func (given function), a (starting point of integral), b (ending point of integral), and n (the number of subintervals). Then, it loops over a range of n subintervals and continuously sums using the trapezoid equation and our given parameters to get the approximate answer of the integral for n subintervals. While the trapezoid equation isn't called in the main body, it is called by another very important definition: the approximator.

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

This function definition is what passes the function, a, and b to the trapezoid_function. It also passes the true value of the integral (calculated through some other method) and the precision to which we want to know the accuracy of the approximation. It tracks the number of subintervals, the approximated sum of that specific subinterval number, and the error from the true answer. The approximator function will also stop once an accuracy to whatever given significant figure is reached. I changed up my method of calculation for this cutoff from notebook 3. I originally used half the final decimal place of wanted degree of accuracy (ie. 10^<sup>-4</sup> accuracy meant 0.00005 was passed). However, this didn't always result in the accuracy we were interested in, so I found a nice slideshow from Illinois.edu that gave me the corrected version we see above. It looks to find when the relative error is less than or equal to 10<sup>-n + 1</sup>, where n is the decimal significant figure. An example of this is if relative error is 10<sup>-2</sup> then the approximation of x has at most three significant figures. The approximator calls on error_calc to find this relative error and then compare it to the passed sig_fig value. Since we were interested in accuracy to the 10<sup>-6</sup> degree, I passed 0.00005 as our sig_fig value to get a value within that error. Until it finds that n, it will double the number of subintervals each loop. The result of this is the following table:

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

As we talked about in the first section, we discovered a better integration method through averaging the leftpoint and rightpoint methods. Now, consider Simpson's rule which fits a parabola using three points from a combination of the midpoint and trapezoid methods. This improved our efficiency, but we can improve it even more by adjusting the sample points themselves. By mapping our original [a, b] range to [-1, 1], we can seek optimal x-coordinates at which to evaluate the given function through a weight function. This range is typically the best for Legendre and Chebyshev polynomials, which we will talk about in the next section. In order to complete this mapping, we will use the following u-substitution:

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

To investigate this idea further, we need to understand Legendre polynomials. Legendre polynomials are a system of complete and orthogonal polynomials with a wide range of mathematical applications. Discovered in 1782 by Adrien-Marie Legendre as the coefficients of the series expansion of Newtonian potentials for electrostatics and gravitation, they are dependent on order and can be represented as P<sub>n</sub> and generated by Rodrigues' formula

```math
P_n(x) = \dfrac{1}{2^nn!}\dfrac{d^n}{dx^n}(x^2-1)^n
```

However, the more important aspect of this set is their orthogonality property. We define this with
```math
\int_{-1}^1 P_m(x)P_n(x)\,dx = \dfrac{2}{2n+1}\delta_{mn}
```
where $\delta_{mn}$ is the Kronecker delta which equals 0 when m &ne n and 1 when m = n. 

This property is what will lead us to using these polynomials for our Gaussian quadrature. However, first let us investigate this property with some code. We were tasked with creating a 4x4 grid of graphs showcasing these polynomials plotted [-1,1] along with their combination of products. To plot the individual polynomials, I just cycled through two for loops and increased the n value of the polynomial for each loop. The legendre was found using the scipy library. I created a separate function definition For the product of the two individual polynomials:

```python
'''
    Definition: Legendre_multiplier
    Parameters: i, j (trackers of loop), x_values (array of values for finding values of y)
    Description: For every loop of the subplots, this calculates the product and returns a list of y
        values to be plotted and the values of the integrated product.
'''
def Legendre_multiplier(i, j, x_values):
    p1 = legendre(i+1)
    p2 = legendre (j+1)

    u_x_values = u(x_values, min(x_values), max(x_values))

    y = p1(u_x_values) * p2(u_x_values) #  Evaluated function for plotting
    integrand = p1 * p2 # Unevaluated function for integrating

    scale = du(min(u_x_values), max(u_x_values)) 

    integrated_value, _ = sp.integrate.quad(integrand, min(u_x_values), max(u_x_values))
    integrated_value = integrated_value * scale 

    return y, integrated_value
```
which returns the plotted values of legendre(k)(x) for the y-axis as well as the integration of the product. Let us check if our previous definition with the Kronecker delta is true and see what these graphs look like. 

![Image](CompProject1Figure.png)

*Fig. 2) This is an image depicting the plots of two Legendre polynomials along with their product, P(i)\*P(j). The subgraphs start at P(1) and P(1) and increase in i and j along the rows and columns to create a 4x4 group of subplots ranging from P(1) to P(4).*

The resulting integration from -1 to 1 for this is given in the following table:

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

This orthogonality property is what will allow us to define the Legendre's roots and weights as points for Gaussian Quadrature. It acts similar to Simpson's rule with how it weights points differently depending on the optimal method of minimizing error. The orthogonality property is important because Legendre polynomials are orthogonal to all polynomials of a lower degree. Using the roots as points means we can minimize this error because it is proportional to our polynomial roots, meaning any polynomial with a lower degree is 0 because of this orthogonality principle. Kind of similar to the step we took to get to the trapezoid rule by minimizing the leftpoint and rightpoint errors. Therefore, we can calculate exact integrations for polynomials for degrees 2n-1 or less. From Dr. Reid's NumericReps file, we can write the approximation for the integral as the following summation:

```math
\int_{-1}^{1} \mathrm{d}x\, f(x) \approx \sum_{i=1}^N c_{N,i} f\left(x_{N,i}\right)
```

where the points $`x_{N,i}`$ are the roots of the N<sub>th</sub> order Legendre polynomial, and the weights are given by the following integral:

```math
c_{i,n}=\frac{1}{P_n^{\prime}(x_{N,i})}\int_{-1}^1\frac{P_n(x)}{x-x_{N,i}} \mathrm{d}x
```

As an additional source, the following diagram depicts the difference in weighting for points between the trapezoid, Simpson's, and Gaussian Quadrature rules:

![Image](CompFig2.png)
*Fig. 3) The interpolation polynomials for various quadrature rules.*

An example of some of these roots and weights are given in the table for the following Legendre polynomials one through four:

        P(x)                                              Roots                                            Weights
    0     1                                              [0.0]                                              [2.0]
    1     2          [-0.5773502691896257, 0.5773502691896257]                                         [1.0, 1.0]
    2     3     [-0.7745966692414834, 0.0, 0.7745966692414834]  [0.5555555555555558, 0.8888888888888883, 0.555...
    3     4  [-0.8611363115940526, -0.3399810435848563, 0.3...  [0.3478548451374538, 0.6521451548625462, 0.652...


### Extension 2: Challenge Boogaloo

I chose to answer the following question:\
Why are the optimal points for an $N$ order Gaussian quadrature the zeros of $P_N$?

My work closely follows the PDF AM205: Gaussian quadrature in my sources (Note: I also included a PDF of my original writing before I turned it into the following markdown). Let us start by defining a generic orthogonal polynomial set (not Legendre yet!) such that the inner product takes the form
```math
&ltp,q&gt = \int_a^b p(x)q(x)w(x)\,dx
```
where w(x) is an arbitrary weight function and for a set of orthogonal polynomials {u<sub>1</sub>, u<sub>2</sub>, etc}:

```math
\langle u_i, u_j \rangle = 0 \quad \text{for } i \neq j
```
Now that we have this paper's definition, we can do the proof.

<ins>Proof.</ins>

Suppose w(x) is an arbitrary weight function and A = {u<sub>1</sub>... u<sub>n</sub>} is the associated orthogonal polynomial set that spans all polynomials of degree less than or equal to n. Let u<sub>n+1</sub> be the associated orthogonal polynomial with degree n+1. Now, consider a monomial, or single term expression, x<sup>L</sup> where L $\leq$ n. Now, since orthogonal polynomials form a basis, we can write x<sup>L</sup> in terms of A such that

```math
x^L = \sum_{i=0}^N \gamma_i u_i(x)
```
where $\gamma_i$ is a scalar.

Then, we can rewrite 
```math
&ltp,q&gt = \int_a^b p(x)q(x)w(x)\,dx
```
as

```math
\implies &ltp,q&gt = \int_a^b x^Lu_{n+1}(x)w(x)\,dx
```

```math
\implies &ltp,q&gt = \int_a^b (\sum_{i=0}^N \gamma_i u_i(x))u_{n+1}(x)w(x)\,dx
```
where we can rewrite the integrand with our defined inner product if we allow p(x) = u<sub>i</sub>(x) and q(x) = u<sub>n+1</sub>(x) such that
```math
\implies \sum_{i=0}^L \gamma_i \langle u_i, u_{n+1} \rangle = 0
```
since u<sub>n+1</sub>(x) is orthogonal to the set A, and every u<sub>i</sub>(x) is an element of that set.

Now, let the points of our quadrature x<sub>0</sub>(x), x<sub>1</sub>(x),...,x<sub>n</sub>(x) be our roots for the polynomial u<sub>n+1</sub>(x)  and define our weights as

```math
w_k = \int_a^b P_k(x)w(x)\,dx
```
where P<sub>k</sub>(x) is the k<sup>th</sup> basis polynomial. Then, we can write an arbitrary polynomial (of degree less than or equal to 2n +1) as the divisible decomposition 
```math
f(x) = p(x)u_{n+1}(x) + r(x)
```
where p and r are polynomials of, at most, degree n. The important term to keep note of is the remainder term r(x). Thus, we can write an integral of f(x) such that

```math
I[f]= \int_a^b f(x)w(x)\,dx
\implies \int_a^b (p(x)u_{n+1}(x) + r(x))w(x)\,dx
= \int_a^b (p(x)u_{n+1}(x)w(x) + \int_a^b r(x))w(x)\,dx
```
where the first term can be written as an inner product and since p(x) is a polynomial of at most degree n, then it is orthogonal to u<sub>n+1</sub>(x), leaving us with
```math
\implies \int_a^b r(x))w(x)\,dx
```
Applying our quadratic algorithm, we find
```math
Q[f] = \sum_{k=0}^n w_kf(x_k)  = \sum_{k=0}^n w_k(p(x_k)u_{n+1}(x_k)+r(x_k))
    = \sum_{k=0}^nw_kr(x_k)
```
These p terms vanish because we let x<sub>k</sub> be the roots of u<sub>n+1</sub>(x), so
```math
u_{n+1}(x_k) = 0 \implies p(x_k)u_{n+1}(x_k) = 0
```
Q[f] approximates the integral of f(x)w(x), but given our f(x) definition we find that
```math
Q[f] = \sum_{k=0}^n w_kf(x_k) = \int_a^b r(x))w(x)\,dx = I[f]
```
Therefore, this is an exact integration for an f $\in$ **P**<sub>2n+1</sub>. This proof is complete once you let w(x) = 1 for Legendre polynomials (or u<sub>k</sub> = P<sub>k</sub>). The PDF notes that you can do a similar method to derive other quadrature polynomial estimates. Therefore, the optimal points for a Guassian quadrature are the zeros of P<sub>n</sub> because it approximates to an exact integration. $\blacksquare$

## Languages, Libraries, Lessons Learned

The primary language for this assignment was Python where we used the libraries scipy, numpy, and pandas. We have consistently worked in Python from the beginning of this module to the final project. The pandas library was especially useful in creating tables and organizing information. I also learned how to create subplots. I knew you could make a 2x2 grid of plots, but I didn't know you could make sizes up to 4x4, so that was neat. Additionally, I enjoyed extension 2 because I am in linear algebra 2 right now, and I saw some similar techniques to what we are doing in that class right now.

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
https://en.wikipedia.org/wiki/Legendre_polynomials#Rodrigues'_formula_and_other_explicit_formulas (Legendre polynomials)
https://faculty.gvsu.edu/boelkinm/Home/ACS/sec-5-6-num-int.html (Left, right, mid figure)

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
