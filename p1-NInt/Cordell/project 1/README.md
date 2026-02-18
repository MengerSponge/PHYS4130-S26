# Introduction

This project demonstrates and explores several numeric integration techniques.

## Background Theory
The first technique explored is the trapezoid rule. The trapezoid rule takes the average of the leftpoint and rightpoint rules. These methods incorporate Riemann sums which split the area under the curve into N rectangles and adds up the area of each rectangle. The height of each rectangle is typically taken at the left, right, or middle of the rectangle--- hence the names of the rules: leftpoint, rightpoint and midpoint. As one might imagine, the estimate becomes more accurate at higher values of N. Specifically, as N increases the error of leftpoint or rightpoint decreases proportionally to $\frac{1}{N}$ and the error of midpoint or trapezoid decreases proportionally to $\frac{1}{N^{2}}$.

Simpson's method is a weighted average of the midpoint and trapezoid rule where midpoint is weighted by a factor of $\frac{2}{3}$ and the trapezoid is weighted by a factor of $\frac{1}{3}$. It is the most accurate out of the methods so far. As N increases, the error reduces by a factor of $\frac{1}{N^{4}}$.

The Gaussian quadrature method is a more refined and much more interesting method than any of the previous methods that revolve around Riemann sums. It utilizes Legendre polynomials to represent the integrand. Recall that Legendre polynomials are an orthonormal set of polynomials as shown in Fig. 1, and so they can be used to represent any polynomial exactly. Furthermore, if the integrand isn't a polynomial, Legendre polynomials can still be utilized to approximate it. The Gaussian quadrature method approximates the integrand using Legendre polynomials and then integrates the constructed polynomial exactly. This is expressed as: 

```math
\int_{-1}^{1} \mathrm{d}x\, f(x) \approx \sum_{i=1}^N c_{N,i} f\left(x_{N,i}\right)                 \qquad(1)
```

where $`x_{N,i}`$ are the sample points chosen at the roots of the Legendre polynomials and the weights are given by

```math
c_{i,n}=\frac{1}{P_n^{\prime}(x_{N,i})}\int_{-1}^1\frac{P_n(x)}{x-x_{N,i}} \mathrm{d}x                  \qquad(2)
```
<br>
<br>
<br>
![Legendre Plot](../Project%201/LegendrePolynomials.png)
Figure 1: Legendre polynomials subplots of $`P_i`$, $`P_j`$, and $`P_i\cdot P_j`$. There are 16 subplots of Legendre polynomials for i and j in the range 1-4. The area under the curve can be observed for each of the i=j subplots on the diagonal as well as the i $`\neq`$ j subplots off diagonal. The areas under the i = j curves look to be approximately 1. Meanwhile, the areas under the i $`\neq`$ j seem to be about 0. This is the visualization of the integral of $`P_i\cdot P_j`$ from -1 to 1 equals the kroenecker delta function.

## Instructions
Code.py is a script that demonstrates the numerical integration techniques discussed above.

To run the code.py script open the command prompt from the directory you have the it in:
```cmd
conda activate [the name of your environment]
python code.py
```
It outputs the estimated result of the integral 

```math
I = \int_0^2 \mathrm{d}x\, \sin^2\left(\sqrt{100x}\right)                   \qquad(3)
```

using the trapezoid rule along with the error and number of subintervals. Then, it outputs the estimated result using the Simpson rule. Lastly, it outputs these results using the Gaussian quadrature.

# Procedure
Here are some highlight-worthy segments of code:

```python
def leftpoint(f,a,b,N):
    mysum = 0
    h = f
    x_array = np.linspace(a, b, N)
    w = (b-a)/(N-1)
    A_array = h(x_array)*w
    mysum = np.sum(A_array)
    mysum = mysum - A_array[N-1]      
    return mysum

def rightpoint(f,a,b,N):
    mysum = 0
    h = f
    x_array = np.linspace(a, b, N)
    w = (b-a)/(N-1)
    A_array = h(x_array)*w
    mysum = np.sum(A_array)
    mysum = mysum - A_array[0]     
    return mysum 

def trapezoid(f, a, b, N):
    return 0.5* (leftpoint(f, a, b, N) + rightpoint(f, a, b, N))
```

The leftpoint and rightpoint functions segment the region being integrated over into N-1 equal width rectangles. The height is calculated using the value of the integrand when x equals the left or right end of each rectangle. The area of each rectangle is calculated and summed to approximate the integral.

```python
def quad(a, b, N):
    roots, weights = sp.special.roots_legendre(N)
    x = ((b-a)*roots/2)+(a+b)/2
    dx_over_du = 2/(b-a)
    return dx_over_du* np.sum(weights*sin(x))
```

The Gaussian quadrature function selects the specific sample points and their relative weights using the scipy special functions library. The integral is parameterized from -1 to 1 because that is where the Legendre polynomials are defined. Then the results are summed over following Eq. (1).

# Analysis
Table 1 compares the approximated results of Integral (3) using the trapezoid rule and the Gaussian quadrature at subintervals incremented by powers of 2. The Gaussian quadrature estimate reaches e-15 error at N=16, whereas the Trapezoid rule increases proportionally to $\frac{1}{N^{2}}$ with e-4 error at N=512. Although over a very very large number of subintervals the trapezoid rule can theoretically reach higher precision, the Gaussian quadrature is much more practical and efficient.

Table 1: Trapezoid and Gaussian quadrature estimates of Integral (3) at N subintervals.

|     N     | Trapezoid Estimate         | Trapezoid Error               | Gaussian Quad Estimate        | Gaussian Quad Error           |
|-----------|----------------------------|-------------------------------|-------------------------------|-------------------------------|
| 2         | 0.9999753123966121         | 0.0057272304291179355         | 0.04681225905124554           | 0.9588902837744845            |
| 4         | 1.421250657401599          | 0.41554811457586904           | 1.4373009028449348            | 0.4315983600192048            |
| 8         | 1.0250832979205577         | 0.01938075509482773           | 1.045246394223079             | 0.03954385139734895           |
| 16        | 0.9467370655877161         | 0.05896547723801393           | 1.0057025428257274            | 2.6645352591003757e-15        |
| 32        | 0.9770282378161866         | 0.028674305009543377          | 1.0057025428257247            | 5.329070518200751e-15         |
| 64        | 0.9976716296671655         | 0.008030913158564523          | 1.0057025428257258            | 4.218847493575595e-15         |
| 128       | 1.003657838280232          | 0.0020447045454980994         | 1.005702542825726             | 3.9968028886505635e-15        |
| 256       | 1.00519113947139           | 0.0005114033543400787         | 1.0057025428257265            | 3.552713678800501e-15         |
| 512       | 1.0055749301202703         | 0.00012761270545968983        | 1.005702542825715             | 1.509903313490213e-14         |

# Extension 1
As shown in Fig. 2, the Gaussian quadrature reaches e-4 error at N=16384, whereas after the substitution it reaches e-15 error at N=16. Furthermore, there is a trend where the error decreases by $\frac{1}{2}$ as N increases as powers of 2 before the substitution. Following this trend, we would expect for it to reach e-10 error around N= $`{2^{31}}`$ . I believe the reason that this is the case has to do with the fact that it is integrating over a singularity, and applying the substitution gets rid of the singularity. Without getting too much into the theory behind Gaussian quadratures (mainly because I don't know all the details lol), this makes sense because polynomials do not have singularities. After all, the Gaussian quadrature approximates the integrand as a polynomial.

![Extension plot](../Project%201/Extension1.png)
Figure 2: Output of Extension 1 code. It shows the number of subintervals, Gaussian quadrature estimate, and error before and after applying the substitution.

# Questions

## Attribution
I used the first two notebooks for the trapezoid rule and Simpson's rule. I was already familiar with these enough that I didn't really have to use any extra resources to learn about them for this project. However, I had never heard of the Gaussian quadrature before. I was able to obtain a fair conceptual understanding of what it does by asking you about it in class.

## Timekeeping
Approximate values 

Sunday: 2 hours
I defined most of the functions I would need for the project and implemented the Trapezoid rule. Then I set up LaTeX, but ended up using markdown in the end.

Monday: 0 hours

Tuesday: 1.25 hours
I spent all of this time learning how to use github properly: pull requests, merging, etc.

Wednesday: 2.5 hours
I implemented the Gaussian quadrature and plotted the Legendre polynomials.

Thursday: 2 hours
I started working on the write up. I explained the background theory for the trapezoid rule and Simpson's rule. I wrote the instructions to run code.py and I picked out and highlighted segments of code.

Monday(2/17): 4 hours
I had an exam in an online class Sunday and another exam today (monday) so I spent the weekend and this morning preparing for those. I wrote the code to perform Extension 1. I talked about the background theory for the Gaussian quadrature in my write up. I explained the previously highlighted segments of code in my write up. I wrote the Analysis, Extension 1, and Questions sections of my write up.

## Languages, Libraries, Lessons Learned
Language: I used python from start to finish.
Libraries: I used numpy, scipy, and pylab. I would say these are all pretty remarkable in their own rights. Numpy and pylab make life easy when doing math and plotting. I use these very frequently. Meanwhile, scipy made it a walk in the park to find the roots and their weights of the Legendre polynomials. I do not use scipy as frequently as numpy and pylab, but I have used it a few times before. Overall, I would say numpy makes life a little bit easier on a frequent occasion, whereas scipy (at least to my experience) makes life a whole lot easier but only once in a while.
Lessons learned: The biggest thing I learned about from this project was the Gaussian quadrature. I think I stated this previously, but I had never heard of it before and now I feel like I have a fair conceptual understanding of what it does, and I wrote code that pulls it off. Another thing I think is worth mentioning, although it isn't super physics related or anything crazy, is that I learned how to code in markdown e.i. formatting tables, importing figures, fence blocks, and niche syntax. 