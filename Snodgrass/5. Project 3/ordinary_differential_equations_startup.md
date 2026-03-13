---
meta:
    author: Austin Reid
    topic: ODE Project
    course: TN Tech PHYS4130
    term: Spring 2026
---

[Symplectic Integrators (libretexts)](https://math.libretexts.org/Bookshelves/Differential_Equations/Numerically_Solving_Ordinary_Differential_Equations_(Brorson)/01%3A_Chapters/1.07%3A_Symplectic_integrators)

# The Assignment

After several class sessions spent on numeric integration, differentiation, and two ODE solving techniques, it's time to do a little more.
In this next project, you will use ODE solvers to explore a few different systems.

Feel free to compare [SciPy's integrators](https://docs.scipy.org/doc/scipy/reference/integrate.html) to your Euler's Method and your RK2.

 1. For a simple harmonic oscillator *not* at rest, what shape will it trace out in phase space?
 3. How does a SHM evolve with a linear damping term?
 4. How does the total mechanical energy (T+V) change over time with different solvers (see below) in these two cases?

## The Deadline

> [!IMPORTANT]
> March 26 (tentatively)

Rough draft (code + outline): March 12

## Program specification:

You will need to compare at least 3 different integrators: One must be symplectic, one must be RK4(5), and the last should be another SciPy integrator.

For the symplectic integrator, you will need to write your own.
The simplest symplectic algorithm is a Verlet integrator.
As an extension you may want to implement a higher order symplectic integrator.

### Extensions (a.k.a. challenge questions)
> [!CAUTION]
>  Do not implement or explore these until you have the base program and report put together.

#### Symplectics Deep Dive: Programming
Generate a method to calculate a phase space volume under time evolution.
The fast way to accomplish this is with four points in (q,p) space and a method to calculate the area of their quadrilateral as a function of time.
This will miss details as you go to longer time, because the shape may stop being rectangular.
Instead, create a small meshgrid of (q,p) points and calculate the area of the concave hull that bounds those points.

 1. How does phase space area (Total Mechanical Energy) evolve in the damped and undamped SHM scenerios?
 2. How do different integrators affect this variation in area?

#### Symplectics Deep Dive: Calculation
Similar to the LibreTexts derivation above, demonstrate that your higher order symplectic integrator conserves phase space.

#### Virial Theorem
Instead of simulating SHM, study one of two multi-particle systems.

 1. Gravity
 2. ball and spring model of a solid

For Gravity, make a function that applies gravitational force between masses.
Write a function that can handle N particles interacting gravitationally.
Start the particles with enough energy to get motion and orbits, but not so much that the collection simply explodes.

For the ball and spring model, write a function that applies spring interactions to an array of masses.
Give the particles enough energy to vibrate, but not so much that they overlap, explode, or disintegrate.

Plot T(t), the kinetic energy of the system over time.
Compare that to V(t), the total potential energy over time.
As N increases, what is the relationship between T(t) and V(t)?
Specifically, what is $`\langle T\rangle`$ and $`\langle V\rangle`$?

## The Writeup

### Address the following questions in your writeup

### Attribution
What resources did you use on this assignment? People, websites, books, etc.

### Timekeeping
How long did you spend on this assignment? If you didn't keep an accurate log, an estimate is fine.

### Languages, Libraries, Lessons Learned
 1. If you changed languages, please remark on that. Otherwise I'll assume you started and finished your work in Python.
 2. What libraries did you use in your submission? Were any of them remarkable? Great to use, super annoying to use, etc?

> [!NOTE]
> This section probably shouldn't more than a few sentences long. Record what you learned and move on!

## Generative AI

> [!IMPORTANT]
> With full attribution and logged prompts, you may use AI on the programming task. You shouldn't, but you may.

> [!WARNING]
> You may _not_ use generative AI on the writeup.

# The Submission

## Code directory
A directory titled "code" that contains your program files and makefile (if used), or a ReadMe.md file that gives sufficiently detailed instructions for me to build and execute your code.
If you implement one or both extensions, make sure it's easy to run the base case as well as the extension you wrote!

## Writeup directory
A directory titled "writeup" containing 
 1. your writeup file: a well-structured markdown file OR a LaTeX document (and its resulting PDF)
 2. Figures (static, PNG or PDF) of the following
 4. Extensions/challenge questions: If you implemented any of the challenges, please include it or them in your report.
 5. *ONLY REQUIRED WITH GENERATIVE AI:* A file titled genAI.md, containing the queries and responses from the AI tool you used. You **MUST** list the following details:
    1. the AI model you used (If you used multiple, indicate which one returned which prompt)
    2. The prompts you supplied
    3. The model's response to the prompt
    4. What part of the generated code you used
    5. What's wrong with the generated code

> [!NOTE]  
> This writeup will need to have figures. It's probably easiest to just add those files to your writeup directory and link them with relative markdown links.

> [!TIP]
> Your figures need to be captioned, which is probably easiest to accomplish with the [FigCaption](https://www.w3schools.com/tags/tag_figcaption.asp) attribute.

> [!TIP]
> Once you have a completed draft, book a consultation with [Tech's Writing Excellence Studio](https://www.tntech.edu/cas/english/west.php)

# Grading
 1. Explanation of undergirding theory
 2. Quality of workflow (including commit history)
 3. Documented implementation of relevant algorithms
 4. Clear and professional presentation of work
