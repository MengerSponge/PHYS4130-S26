---
meta:
    author: Austin Reid
    topic: Numeric Integration Project
    course: TN Tech PHYS4130
    term: Spring 2026
---

# Numeric Integration

So far you've worked through two worksheets from NCSU's computational physics course and a problem set from a Berkeley Computational Methods textbook.
You're ready for a little more sophistication.

# The Assignment

## Generative AI

> [!IMPORTANT]
> With full attribution and logged prompts, you may use AI on the programming task. You shouldn't, but you may.

> [!WARNING]
> You may _not_ use generative AI on the writeup.

## The Deadline (submissions are via pull requests to the relevant directory)

### Rough Draft (primary assignment and writeup, no extensions)

> [!IMPORTANT]
> Beginning of class on ???

### Final Submission

> [!IMPORTANT]
> Beginning of class on ???

## Programming specification:

### Adaptive integration: Trapezoidal Rule

Consider the following integral:

```math
I = \int_0^2 \mathrm{d}x\, \sin^2\left(\sqrt{100x}\right)
```

It's easy to implement a trapezoidal rule approximation, and with a high enough $N$ you can generate an arbitrarily accurate solution.
It's wasteful to just crank up $N$ to some huge value though.
Instead, it's helpful to consider "how good" the approximation is at any particular point.

> [!IMPORTANT]
> Calculate the integral to an approximate accuracy of $`\epsilon=10^{-6}`$, which means it's correct for 6 digits after the decimal point.

Start with one single integration slice and work up from there to two, four, eight, and so forth.
For each value of the number of slices $`N`$: your program should print out the number of slices, its estimate of the integral, and its estimate of the error on the integral.

> [!TIP]
> The value is $`\approx 1.0057`$

### Gaussian quadrature

You saw the performance gains made possible by moving from the trapezoidal rule to Simpson's rule.
Simpson's rule effectively fits a parabola to each set of three points, and uses that to approximate the integral.
This is a more sophisticated *quadrature*. (mostly archaic, but still used in the following way)

Simpson's rule reduces error by cleverly weighting the contributions of points that you sample along the interval you're integrating.
Without a huge amount of justification, it seems reasonable that you can get a better approximation by adjusting your sample points themselves.

#### Standardized limits of integration

It turns out it's a bit easier to only consider integrals in the interval $`[-1,1]`$ 

> [!NOTE]
> It's possible to extend these results to other intervals, but that takes us further from *neat* and much closer to *tedious*.
> The intervals don't even have to be finite! $`[0,\infty)`$, for example, is manageable.)

You generally won't be lucky enough to attempt an integral over this interval.
Instead, you're more likely to see an integral from $`[a,b]`$.
You only need to map $`[a,b]`$ to $`[-1,1]`$, a simple task for the right u-substitution.
Let

```math
u=\frac{2x-a-b}{b-a}.
```

> [!IMPORTANT]
>  1. Verify that u(x) maps $`[a,b]`$ to $`[-1,1]`$.
>  2. Calculate $`\mathrm{d}u`$ and make sure you remember to use it in your program!

#### Legendre polynomials

You've seen legendre polynomials in Quantum Mechanics.
They're "orthogonal polynomials", but a brief refresher will be helpful.

Make a figure with a 4x4 grid showing the following subplots: $`P_i`$, $`P_j`$, and $`P_i\cdot P_j`$

|       | $`P_1`$                  |   $`P_2`$                   |  $`P_3`$                    |  $`P_4`$  |
| :---: | :----:                 | :----:                 | :---:                  | :---:                  |
| $`P_1`$    | $`P_1`$, $`P_1`$, $`P_1\cdot P_1`$ | $`P_1`$, $`P_2`$, $`P_1\cdot P_2`$ | $`P_1`$, $`P_3`$, $`P_1\cdot P_3`$ | $`P_1`$, $`P_4`$, $`P_1\cdot P_4`$ |
| $`P_2`$    | $`P_2`$, $`P_1`$, $`P_2\cdot P_1`$ | $`P_2`$, $`P_2`$, $`P_2\cdot P_2`$ | $`P_2`$, $`P_3`$, $`P_2\cdot P_3`$ | $`P_2`$, $`P_4`$, $`P_2\cdot P_4`$ |
| $`P_3`$    | $`P_3`$, $`P_1`$, $`P_3\cdot P_1`$ | $`P_3`$, $`P_2`$, $`P_3\cdot P_2`$ | $`P_3`$, $`P_3`$, $`P_3\cdot P_3`$ | $`P_3`$, $`P_4`$, $`P_3\cdot P_4`$ |
| $`P_4`$    | $`P_4`$, $`P_1`$, $`P_4\cdot P_1`$ | $`P_4`$, $`P_2`$, $`P_4\cdot P_2`$ | $`P_4`$, $`P_3`$, $`P_4\cdot P_3`$ | $`P_4`$, $`P_4`$, $`P_4\cdot P_4`$ |

Look at the $`P_i\cdot P_j`$ plots, and convince yourself that if $`i=j`$, $`\int_{-1}^{1} P_i\cdot P_j \mathrm{d}x=\delta_{i,j}`$.
The Kroneker delta function $`\delta_{i,j}`$ is 1 when $`i=j`$ and 0 when $`i\neq j`$.

For deep numerical analysis reasons that you're welcome to read up on (and are encouraged to address in your report if you like), the optimal points for an $N$ order Gaussian quadrature are the zeros of that order of Legendre polynomial.
You won't need to define your own Legendre polynomial generating function: a standard numerical library implementation is fine.

#### The algorithm

```math
\int_{-1}^{1} \mathrm{d}x\, f(x) \approx \sum_{i=1}^N c_{N,i} f\left(x_{N,i}\right)
```

Where the points $`x_{N,i}`$ are the roots of the Nth order legendre polynomial, and the weights are given by the following integral:
```math
c_{i,n}=\frac{1}{P_n^{\prime}(x_{N,i})}\int_{-1}^1\frac{P_n(x)}{x-x_{N,i}} \mathrm{d}x
```

Note that the sample points and the relative weights are *purely* functions of the Legendre polynomials!
Gauss calculated several of them, and their values haven't changed since.
It might be a huge pain to calculate them, but you only have to do it *one time, EVER*.

For now, you can use one of the standard libraries to get the roots of $`P_N`$:
```python
import scipy as sp
roots, weights = sp.special.roots_legendre(N)
```

Conveniently, scipy's roots_legendre function returns the weights too.

### Extensions (a.k.a. challenge questions)

> [!CAUTION]
> Do not implement or explore these until you have the base program and report put together.

Consider the integral

```math
\int_0^2\mathrm{d}y\, \frac{y^2}{\sqrt{2-y}}
```

At the upper limit of integration, the integrand becomes infinite. This integral still converges ($`\sqrt{8192}/15`$), but we have to be very careful with our numeric tools around those pesky infinities.
Gaussian quadrature points are not evaluated at the limits of integration.

 - How many points do you need in your Gaussian quadrature to achieve 10 digits of precision?
 - Apply the change of variable $`y=2\sin^2\theta`$. 
 - How many Simpson's rule points do you need to calculate this to 10 significant figures?
 - Now use Gaussian quadrature to calculate the integral. How many points do you need to achieve that same precision?


### Extension 2: Challenge Boogaloo

This is a numeric analysis/theory question, not a scientific programming task.
Why are the optimal points for an $N$ order Gaussian quadrature the zeros of $P_N$?

<!-- > [!NOTE]
> Convergence, stability, and error will feature heavily in this course.
>  - Will the algorithm described above always terminate?
>  -  -->

## The Writeup

### Address the following questions in your writeup

### Attribution

What resources did you use on this assignment? People, websites, books, etc.

### Timekeeping

How long did you spend on this assignment? If you didn't keep an accurate log, an estimate is fine.

### Languages, Libraries, Lessons Learned

 1. What language did you use for your submission? Is it the same one you started using? If not, why'd you change?
 2. What libraries did you use in your submission? Were any of them remarkable? Great to use, super annoying to use, etc?

> [!NOTE]
> This section probably shouldn't more than a few sentences long. Record what you learned and move on!

# The Submission

Your submission should explain what your code does, with highlighted snippets as appropriate.
It should address the broader questions asked above, and some of your specific findings too.
When in doubt, include it!

Remember the vocabulary you've learned regarding efficiency, precision, and accuracy.
Speed isn't exactly the same thing as efficiency, but they're closely related.
Consider using timing information as you evaluate your implementation of these different types of quadrature.

> [!IMPORTANT]
> Make sure you run spell check on your text file.

## Code directory

A directory titled "code" that contains your program files and makefile (if used), or a ReadMe.md file that gives sufficiently detailed instructions for me to build and execute your code.
If you implement one or both extensions, make sure it's easy to run the base case as well as the extension you wrote!

## Writeup directory

A directory titled "writeup" containing 
 1. your writeup file: a well-structured markdown file OR a LaTeX document (and its resulting PDF)
 3. Figures (static, PNG or PDF)
 4. Extensions/challenge questions: If you implemented any of the challenges, please include it/them in your report.
 5. *ONLY REQUIRED WITH GENERATIVE AI:* A file titled genAI.md, containing the queries and responses from the AI tool you used. You **MUST** list the following details:
    1. the AI model you used (If you used multiple, indicate which one returned which prompt)
    2. The prompts you supplied
    3. The model's response to the prompt
    4. What part of the generated code you used
    5. What's wrong with the generated code

> [!IMPORTANT]
> Your figures must render as part of your report.

> [!TIP]
> Your figures need to be captioned, which is probably easiest to accomplish with the [FigCaption](https://www.w3schools.com/tags/tag_figcaption.asp) attribute.

# Grading
 1. Explanation of undergirding theory
 2. Quality of workflow (including commit history)
 3. Documented implementation of relevant algorithms
 4. Clear and professional presentation of work

# Acknowledgements
This assignment is begged, borrowed, and stolen from the following sources:
 -  [Prof. Mark Newman's course at UMich](https://websites.umich.edu/~mejn/courses/2012/phys411/homeworks/homework4.pdf)
 -  [Jeffrey Tatum's Numerical Methods LibreText](https://phys.libretexts.org/Bookshelves/Astronomy__Cosmology/Celestial_Mechanics_(Tatum)/01%3A_Numerical_Methods/1.15%3A_Gaussian_Quadrature_-_the_Algorithm)
