---
meta:
    author: Emma Krebs
    topic: Numeric Integration Project
    course: TN Tech PHYS4130
    term: Spring 2026
---

# Numeric Integration

## Languages, Libraries, Lessons Learned

The primary language for this assignment was Python where we used the libraries scipy, numpy, and pandas. We have consistenly worked in Python from the beginning of this module to the final project. The pandas library was especially useful in creating tables and organizing information. I also learned how to create subplots. I knew you could make a 2x2 grid of plots, but I didn't know you could make sizes up to 4x4, so that was neat. 

## Timekeeping

As of 2/17/26, I have spent 15 hours or so on this project. 

## Summary of Code

### Trapezoid section

There are several methods when approximating a definite integral over a limit through Riemann Sums. This is when rectangles are formed by dividing the integral [a, b] into N subintervals from [x<sub>0</sub>, x<sub>1</sub>] to [x<sub>N-1</sub>, x<sub>N</sub>], where x<sub>0</sub> is a and x<sub>N</sub> is b. This approximates the area under the curve for a function f(x) when mulitplied by the height (or point on f(x)). For increasing subintervals N our approximation tends to become more accurate (Excluding cases with functions with alternating peaks above and below the x-axis since it can widely change between the number of subintervals depending on the points chosen. However, it will eventually approach the correct approximation). 

Although all these methods will get you to your correct approximation, not all methods are equivalent in terms of efficency. This will effect both our timing, making more complex integration's approximations longer when trying to reach a certain precision or error, and computer resources for each loop run.

There are five main methods that we have discussed so far: leftpoint, rightpoint, midpoint methods, trapezoid rule, and Simpson's rule. Leftpoint, rightpoint, and midpoint methods inherently follow their namesake such that their rectangles align with the left, right, and middle points, respectively, of the top edge of the rectangle. An example of the three are shown in the figure below. The error of the leftpoint and rightpoint methods decrease at a similar rate and are the worst efficency approximations out of the five with the error decreasing at a linear rate. The midpoint method and trapezoid rule are slightly better decreasing at a quadratic proportionality. The best of the five is the Simpson's rule, which combines the weighted sums of the midpoint method and trapezoid rule to get an error proportionality with respect to the fourth power. However, for this code we worked on coding a trapezoid method to solve the following integral:

```math
I = \int_0^2 \mathrm{d}x\, \sin^2\left(\sqrt{100x}\right)
```

To do this, let us go a little more in depth with the trapezoid rule. We mentioned before that the proportionality of decreasing error for the leftpoint and rightpoint methods are approximately the same. Thus, the point of evaluation should not make much of a difference assuming our width of a subinterval is small and since the errors should be roughly the same in magnitude, but opposite in sign, we can make have a better method of integration by averaging both methods. This is where we get our trapezoid rule. That is, to approximate the integral
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

-Talk about use of gaussian quadrature (honeslty not sure whhere to use this two functions as of right now...)
-Talk about how it maps to -1 to 1 and why we care
-You will probably add more to this after doing the challenege section.

### Legendre Polynomials
-Talk about use of 





### Extension

-Challenge problem 1 (talked about 2 in guassian quadrature above)

## Sources

People Used:

Cricket/Cordell recommending how to make 4x4 subplot table look better. Dr. Reid is also another contributor after I stole your code from the project's markdown file to make my header look nice.

Websites Used:

https://stackoverflow.com/questions/455612/limiting-floats-to-two-decimal-points (For round function)
https://www.desmos.com/calculator/d9rmt4wfoa (Checked trap rule against)
https://www.geeksforgeeks.org/python/printing-lists-as-tabular-data-in-python/ (Table building)
https://cs357.cs.illinois.edu/textbook/assets/slides/03-Errors.pdf (New method for error relation to sigfig)
https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.legendre.html (Legendre polynomials)
https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.set_xlabel.html (More on axes and matplotlib)
https://personal.math.ubc.ca/~israel/m210/lesson19.pdf (Double checked my error for the integration methods)

Books Used:

Python All-In-One for dummies by John C. Shovic, PHD and Alan Simpson. (Mostly emotional support...)
