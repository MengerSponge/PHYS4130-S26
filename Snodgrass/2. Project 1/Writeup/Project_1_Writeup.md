---
meta:
    Author: Alec Snodgrass
    Topic: Numeric Integration Project
    Course: TN Tech PHYS4130
    Term: Spring 2026
---
# Numeric Integration Project

## Trapezoid approximation

The trapezoidal approximation was developed in earlier notebooks and builds off the left endpoint, right endpoint, and midpoint rules. Those three rules were studied regarding Riemann sums, a building block for integration. Using the Riemann sum method, integrals can be approximated by taking finer and finer slices of the area under the function's plot. 

Where the left endpoint, right endpoint, and midpoint rules differ from the trapezoid approximation is the form of the 'top' of the 'slice'. The trapezoid uses both the left and the right sides to approximate the same slope as the function. This method cuts out a large portion of over- or underestimations made by square rectangles. 

![Integration trapezoids](https://upload.wikimedia.org/wikipedia/commons/a/a8/Calkowanie_numeryczne-metoda_trapezow.png)

The code below implements the function:
```math
  I_T = \frac{h}{2}\left[ f(a) + f(b)\right] + \sum_{i=1}^{N-1} f(a+ih) h
```

```python
def traprule(f, a, b, N):                   
    dx = (b-a) / N
    l = np.linspace(a, b-dx, N)           # This calculates the left endpoint x-value
    r = np.linspace(a+dx, b, N)           # This calculates the right endpoint x-value
    return np.sum( (f(l)+f(r)) * dx ) / 2 # Averages left and right rules to form a trapezoid
```
This function divides the range into N subsections, finds the left and right endpoints for each subsection, finds the value at both endpoints for each subsection, and then calculates the area for each subsection. The area is calculated by averaging the subsection area calculated by the left and right endpoints. The N areas are added and returned as the integral.

### Performance of the traprule
This function was used with a *very* high **N** value to get a good approximation of the integral to compare against later. The integral was approximately:
```math
\int_0^2 \mathrm{d}x\, \sin^2\left(\sqrt{100x}\right) = 1.0057025428257
```
  **Table 1.** Performance of trapezoidal approximation.
  | N         | Integral              | Error                |
  |-----------|-----------------------|----------------------|
  | 2         | 0.79594662529161      | 2.0975591753378 e-01 |
  | 4         | 0.69837008702424      | 3.0733245580114 e-01 |
  | 8         | 1.03497028020385      | 2.9267737378463 e-02 |
  | 16        | 0.94670012038502      | 5.9002422440363 e-02 |
  | 32        | 0.97846523867870      | 2.7237304146686 e-02 |
  | 64        | 0.99790966934232      | 7.7928734830655 e-03 |
  | 128       | 1.00368931558769      | 2.0132272376980 e-03 |
  | 256       | 1.00519511618492      | 5.0742664047077 e-04 |
  | 512       | 1.00557542778294      | 1.2711504245244 e-04 |
  | 1024      | 1.00567074790213      | 3.1794923253337 e-05 |
  | 2048      | 1.00569459308443      | 7.9497409624096 e-06 |
  | 4096      | 1.00570055532724      | 1.9874981440626 e-06 |
  | 8192      | 1.00570204594715      | 4.9687823322486 e-07 |
  | 16384     | 1.00570241860583      | 1.2421955486452 e-07 |
  | 32768     | 1.00570251177073      | 3.1054654403562 e-08 |
  | 65536     | 1.00570253506197      | 7.7634143558214 e-09 |
  | 131072    | 1.00570254088478      | 1.9406036777525 e-09 |
  | 262144    | 1.00570254234049      | 4.8490100823528 e-10 |
  | 524288    | 1.00570254270441      | 1.2097522983368 e-10 |
  | 1048576   | 1.00570254279539      | 2.9993785233273 e-11 |
  | 2097152   | 1.00570254281814      | 7.2484240831727 e-12 |
  | 4194304   | 1.00570254282383      | 1.5623058402525 e-12 |
  | 8388608   | 1.00570254282525      | 1.4077627952247 e-13 |
  | 16777216  | 1.00570254282560      | 2.1493917756743 e-13 |
  | 100000000 | 1.00570254282572      | Control              |

### Significance of these calculations
It is apparent from the table above that the trapezoidal rule is inefficient for approximating complicated integrals. The following sections of the writeup demonstrate and somewhat explain a much more efficient integral approximation technique developed by Gauss. The resulting algorithm provides a efficient, accurate, and fast method for approximating problematic integrals. 



## Legendre polynomial orthogonality

The Legendre polynomials are ubiquitously important as one of the fundamental sets of orthogonal polynomials. The orthogonality allows us to represent functions as a sum of legendre polynomials using an approach resembling the Fourier series expansion, but this time we use the Legendre polynomials as the orthogonal basis functions. 

Before discussing the Gaussain quadrature method in detail, it is important to refresh yourself on the Legendre polynomials and their properties. Any task involving Legendre polynomials will always be made easier if the requirement to calculate them is waived. Thankfully, we were allowed to use a library that handed us the Legendre polynomials straight away.

```python
from scipy.special import legendre
```
This function streamlined the process of calculating polynomials and multiplying them together. Below is a grid of low-order Legendre polynomials. Upon careful inspection of the graphs, it is clear that $`\int_{-1}^{1} P_i\cdot P_j \mathrm{d}x=\delta_{i,j}`$. Meaning that the Legendre polynomials are, in fact, orthogonal. 

<figure>
  <img src="Figures/Legendre_Grid.png" style="width:100%">
  <figcaption><strong>Figure 1.</strong> 4×4 grid of Legendre polynomials and their products.</figcaption>
</figure>

As you may have noticed, the interval for these functions is [-1, 1]. This is an important characteristic that will affect how these polynomials can be used. 

---

## Gaussian quadrature

### Legendre basis and integral setup
Gaussian quadrature uses the Legendre polynomials to construct a series of orthogonal polynomials each with a weighted contribution. The first step when using this method is to define the integrand over the same interval, [-1, 1], that the Legendre polynomials form a basis on. However, most integrals will not be over these bounds, and instead over the bounds [a, b]. Therefore, a simple function conversion is required to use the algorithm. Below, I map the integrand and the bounds to an equivalent function over [-1, 1]. 
```math
u(x) = \frac{2x-a-b}{b-a}
```
Therefore, if x = a:
```math
u(a) = \frac{2a-a-b}{b-a} = -1
```
And, if x = b:
```math
u(b) = \frac{2b-b-a}{b-a} = 1
```
To calculate the integral, we need to account for our change of variables.
```math
u=\frac{2x-a-b}{b-a}
```
```math
u(b-a) = (2x-a-b)
```
```math
u(b-a)+a+b = 2x
```
```math
x = \frac{u(b-a)}{2} + \frac{a+b}{2}
```
```math
\mathrm{d}x = \mathrm{d}u * \frac{b-a}{2}
```
For this problem, [a, b] = [0,2] and therefore 
```math
\mathrm{d}x = \mathrm{d}u \\
x = u+1
```
The function f(x), which is being integrated, is:
```math
f(x) = \sin^2\left(\sqrt{100x}\right)
```
This function is therefore transformed to g(u) so the algorithm can be easily applied
```math
g(u) = \sin^2\left(\sqrt{100(u+1)}\right)
```
The final integral transformation is:
```math
\boxed{\int_0^2 \sin^2\left(\sqrt{100x}\right)\,\mathrm{d}x
\;\Longrightarrow\;
\int_{-1}^1 \sin^2\left(\sqrt{100(u+1)}\right)\,\mathrm{d}u}
```

---
### Algorithm explanation
To approximate the integral, Gaussian quadrature is laid out in the following algorithm
```math
\int_{-1}^{1} \mathrm{d}x\, g(x) \approx \sum_{i=1}^N c_{N,i} g\left(x_{N,i}\right)
```

Where the points $`x_{N,i}`$ are the roots of the Nth order Legendre polynomial, which are given by the **scipy** function  discussed below. The weights are given by the following integral:

```math
c_{i,n}=\frac{1}{P_n^{\prime}(x_{N,i})}\int_{-1}^1\frac{P_n(x)}{x-x_{N,i}} \mathrm{d}x
```

Scipy, arguably, provides a function even more useful than the Legendre polynomials library. 
```python
import scipy as sp
roots, weights = sp.special.roots_legendre(N)
```
This special roots function not only returned the roots of the Legendre polynomial of order N, but also returned the weighting coefficient. Both of which are needed in the summation given above. All that is left to do is multiply and add. 
```python
def algorithm(N):
    roots, weights = sp.special.roots_legendre(N)
    return np.sum(weights * g(roots))
```

---
### Algorithm performance
Using this algorithm, the integral is approximated to an accuracy of $\epsilon=10^{-14}$ with N = 16. The efficiency is remarkable compared to the trapezoidal approximation, which required tens of millions of sub-intervals to get an accuracy of even $10^{-12}$. See the table below for the output of the Gaussian quadrature.

  **Table 2.** Performance of Gaussian quadrature algorithm.
  | N   | Approximation       | Signed Error          |
  |-----|---------------------|-----------------------|
  | 1   | 0.591917938186608   | 0.413784604639128     |
  | 2   | 0.046812259051246   | 0.958890283774491     |
  | 3   | 1.078855067763076   | -0.07315252493733     |
  | 4   | 1.437300902844935   | -0.43159836001919     |
  | 8   | 1.045246394223079   | -0.03954385139734     |
  | 16  | 1.005702542825727   | 9.32587340685131 e-15 |
  | 32  | 1.005702542825725   | 1.19904086659516 e-14 |
  | 64  | 1.005702542825726   | 1.08801856413265 e-14 |
  | 128 | 1.005702542825726   | 1.06581410364015 e-14 |
  | 256 | 1.0057025428257265  | 1.02140518265514 e-14 |
  |2048 | 1.0057025428257322  | 4.44089209850063 e-15 |



## Extension 1
###   Problem definition
The integral given for the challenge is:
```math
\mathrm{Integral\,\, 1:}\,\,\,\int_0^2 \frac{y^2}{\sqrt{2-y}} \, \mathrm{d}y
```
Which converges to
```math
\frac{\sqrt{8192}}{15} \Rightarrow 6.033977866125206
```

Three methods of approximation are used and compared. A table with the results is show at the end. Each method is used to approximate the integral out to 10 digits of precision and their efficiency is compared in the table. 

---
### Change of variables
To make this integral manageable, a change of variables in made below using:
```math
y = 2\sin^2\theta \\
\mathrm{d}y = 4\sin\theta\cos\theta \mathrm{d}\theta
```
**Performing the change of variables**

Numerator:
```math
y^2 = (2\sin^2\theta)^2 = 4\sin^4\theta

```
Denominator: 
```math
2 - y = 2 - 2\sin^2\theta = 2\cos^2\theta \\
\sqrt{2-y} = \sqrt{2}\cos\theta
```

#### **Changing the limits**
When y = 0:
```math
0 = 2\sin^2\theta\\
\theta = 0
```
When y = 2:
```math
2 = 2\sin^2\theta\\
\theta = \frac{\pi}{2}
```
So the new limits are
```math
\theta \in [0, \tfrac{\pi}{2}]
```

#### **Substituting**
The integral becomes
```math
\int_0^{\pi/2}\, \frac{4\sin^4\theta}{\sqrt{2}\cos\theta}\, 4\sin\theta\cos\theta\, \mathrm{d}\theta
```
Which simplifies to:
```math
\mathrm{Integral\,\, 2:}\,\,\,\boxed{  8\sqrt{2}\, \int_0^{\pi/2}  \sin^5\theta\, \mathrm{d}\theta }
```

---
### Simpson's Rule
After the change of variables is made, only a short section of code in necessary to compute Simpson's rule for approximating integrals. (Simpson's rule was developed in a previous notebook. The function is given below but see previous notebooks for the theory behind it.)
```python
def Simpson(integrand, a, b, N):            # Integral on the bounds [a, b]
    I_m = midpoint(integrand, a, b, N)      # Midpoint-rule approx
    I_t = traprule(integrand, a, b, N)      # Trapezoid-rule approx
    return (2*I_m + I_t) / 3 
```
With **N** equal to 25, the approximation is accurate out to 10 digits of precision.

---
### Gaussian Quadrature (for 1 and 2)
The first step in Gaussian quadrature is to map the function onto [-1, 1]. This will need to be done for both the original integral and the change-of-variables integral **(*Integral 1* and *Integral 2*)**. The general change of variables is:

```math
x = \frac{b-a}{2}u + \frac{a+b}{2}
```
```math
dx = \frac{b-a}{2} \mathrm{d}u
```

Therefore, the before-and-after relationship is:
```math
\int_a^b f(x)\,dx = \frac{b-a}{2}\, \int_{-1}^{1}\, f\!\left(\frac{b-a}{2}u + \frac{a+b}{2}\right) \mathrm{d}u
```


**Original Integral**
```math
\int_0^2 \frac{y^2}{\sqrt{2-y}} \, \mathrm{d}y
```
Since the original range is [0, 2], the function is just offset by 1 from our needed range. Therefore, the scaling is a factor of 1. 

```math
y = u + 1
```
```math
dy = dt
```
```math
\frac{b-a}{2} = 1, \quad \frac{a+b}{2} = 1
```

Therefore,
```math
\int_0^2 f(y)\,dy  =  \int_{-1}^{1} f(u+1)\,\mathrm{d}u
```
Which is:
```math
\int_0^2 \frac{y^2}{\sqrt{2-y}}\,\mathrm{d}y  = 
\boxed{ \int_{-1}^{1}\, \frac{(u+1)^2}{\sqrt{1-u}}\,\mathrm{d}u }
```
Briefly, the Gaussian quadrature algorithm is being applied to Integral 1, which has been transformed to fit the range of [-1, 1] requirement. The result of that computation is show in the table below.

---

**Transformed Integral**

After substitution, we obtained

```math
8\sqrt{2}\, \int_0^{\pi/2}\, \sin^5\theta \, \mathrm{d}\theta
```
With these bounds, the relation becomes

```math
\theta = \frac{\pi}{4}(u+1)
```
```math
d\theta = \frac{\pi}{4} du
```
since
```math
\frac{b-a}{2} = \frac{\pi}{4},\,\ \frac{a+b}{2} = \frac{\pi}{4}
```

Therefore,

```math
\int_0^{\pi/2} g(\theta)\,d\theta  =  
\frac{\pi}{4}\, \int_{-1}^{1}\, g\!\left(\frac{\pi}{4}(u+1)\right) \mathrm{d}u
```
Which is
```math
8\sqrt{2}\int_0^{\pi/2} \sin^5\theta\,\mathrm{d}\theta  =
\boxed{ \pi\sqrt{2}\, \int_{-1}^{1}\, \sin^5\!\left(\frac{\pi}{4}(u+1)\right) \mathrm{d}u }
```
Briefly, the Gaussian quadrature algorithm is being applied to Integral 2, which has been transformed to fit the range of [-1, 1] requirement. The result of that computation is show in the table below.

---
### N, Runtime, and Accuracy
  **Table 3.** Comparison of numerical integration methods.
  | Method              | N       | Approximation     | Error                   |
  |---------------------|---------|-------------------|-------------------------|
  | $\frac{\sqrt{8192}}{15}$ | -  | 6.033977866125206 | Actual Value            |
  | Gauss Quadrature 1  | 10,000  | 6.033485338530997 | 0.0004925275942095908   |
  | Simpson's Rule      | 46      | 6.033977866147473 | -2.226663298188214 e-11 |
  | Gauss Quadrature 2  | 10      | 6.033977866125502 | -2.957634137601417 e-13 |


The Gaussian quadrature for the original integral was poor. The algorithm took a very long time (almost four minutes) to run with N = 100,000, and the result was only accurate out to four decimal places! Terribly inefficient. 

The Gaussian quadrature for the change of variables integral was **much** better. The result compiled in 0 seconds and only required N = 9 for an accuracy of 10 significant figures. This beat the Simpson's rule for the changed variable integral, which required N = 25 for the same accuracy. 

## Attribution

- People
  - Dr. Reid answered several clarifying questions throughout the project. Mostly, questions about his intent for our work and the form of our submissions. 
- Websites
  - A few different websites were referenced for coding syntax and function definitions
- Books
  - Who has books anymore?!
- AI
  - Much more modern approach. I asked ChatGPT on several occasions to suggest Python functions and to explain their results, such as condensing plots, formatting graphs, and loop syntax...

## Timekeeping

After the full write-up and extensions are finished, the total time spent will be close to two full day's work. The extension added a full day basically. The amount of writing, typing, and formatting it took was painstaking. 

## Languages, Libraries, Lessons Learned

### Languages
- English mainly, but also Python. I would use C++ or C, but I don't have much experience with Python, and I think it would be good to get some exposure to it. At some point, though, I would like to do it in C++ to refine my proficiency
- I started with Python and never thought about switching.

### Libraries
- I used the standard libraries for math-type coding: 
  - NumPy 
  - SciPy
- I also had to use:
  - matplotlib.pyplot
  - SciPy. Special for the Legendre functions
- There was a REMARKABLE library. I pointed this out in the code as well, but the sp.specail.roots_legendre(N) function was AWESOME! It practically completed the most intense part of the algorithm for me. 
- The pyplot library was difficult to use because I have not had much experience with it, and I don't know the syntax. Typically, that is not a big deal, but for a plotting function, there are a lot of steps and customizations that go into it. 

### Lessons learned
- I got more familiar with the syntax of Python, the libraries and their uses, and some functions for math, plotting, etc. 
- Obviously, I learned about Gaussian quadrature and what a wonderful method of integration it is!
- More interestingly than the program-specific topics, I had to practice program management, timekeeping, version control software, git-specific nuances, and more. I think that this project was a great learning curve for the rest of the course. 