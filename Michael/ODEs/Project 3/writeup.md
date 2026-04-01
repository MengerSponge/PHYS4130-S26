# Project 3
## Introduction
Few areas of math enjoy such a privelleged position in physics as differential equations do. However, for all of their utility in modeling physical system, the average differential equation(s) that we encounter outisde of an introductory course lack any kind of analytic solution. This necessitates the development of numerical methods that can efficiently and accuratley approximate solutions to initial value problems. In this porject, we will examine different algorithms for solving ODEs and see how they solve the simple harmonic oscillator along with some other special cases.

## Algorithms and Theory
We will first examine the different numerical methods that will be used. These functions were written to be highly general in their use case. The functions are able to take in lists of initial conditions and solve a system of coupled first order ODEs. The outputs are structured as lists of solution. For example, if you pass in lists of length n, then the functions output a list with the array of times as its first element, and a sequence of two element lists containing the phase space trajectories for the lists of initial conditions. They are all capable of solving a system of the form

```math
\dot{x} = f(t,x,y)
```
```math
\dot{y} = g(t,x,y)
```
The implementation can roughly be structured as 
```python
X_initial = [1, 2, 3, 4]
Y_initial = [5, 6, 7, 8]

tmin = 0
tmax = 1
nts = 100

def du_dt(t, U) #U = [X, Y] 
    X, Y = U
    dx_dt = f(t, X, Y)
    dy_dt = g(t, X, Y)
    return [dx_dt, dy_dt] #MUST be a list

sol = Method(X_initial, Y_initial, tmin, tmax, du_dt) #generic method

#How to extract information from sol
t = sol[0]

X1 = sol[1][0]
Y1 = sol[1][1]

X2 = sol[2][0]
Y2 = sol[2][1]

X3 = sol[3][0]
Y3 = sol[3][1]

```
With a basic explanation how the methods will all be structured done, we can look at the different method that will be used. The first one is RK4(5). It belongs to a family of solvers known as Runge-Kutta (RK) methods. The explicit derivation is not of interest here, but a brief explanation of how this family of solvers works is useful. They are what's known as prediction corrector methods. Unlike bad solvers such as Euler's method, "predictor-corrector methods improve the approximation accuracy by querying the 𝐹 function several times at different locations (predictions), and then using a weighted average of the results (corrections) to update the state." RK4(5) is unique in that it is an adaptive method. It achieves by chagning the step size through comparing a 4th order step and a 5th order step. A python implementation is simple using the solve_ivp function from scipy.integrate. 

```python
from scipy.integrate import solve_ivp

def RK45(X0, Y0, tmin, tmax, nts, du_dt):
    t_span = (tmin,tmax)
    t = np.linspace(tmin,tmax,nts,endpoint = False)

    #Convert the initial condtions from lists to arrays
    X0 = np.asarray(X0)
    Y0 = np.asarray(Y0)

    #Packaging the solutions for the way I like them. It loops over the initial conditions and       #solves. This does mean that [X1, Y1], [X2, Y2], ... can't be coupled.
    solutions = [t]
    for i in range(len(X0)):
        sol = solve_ivp(du_dt, t_span, [X0[i], Y0[i]], t_eval = t, method = 'RK45')
        solutions.append([sol.y[0], sol.y[1]])

    return solutions
```

The next method used is known as LSODA. It is another adaptive numerical solver for systems of ordinary differential equations. It was developed as part of the ODEPACK library and is designed to efficiently handle both stiff and non-stiff problems without requiring the user to decide which type of solver to use. It automatically switches between two classes of multistep methods depending on the behavior of the system. When the problem appears non-stiff, LSODA uses variable-order Adams predictor–corrector methods, which are explicit multistep schemes that are efficient for smooth solutions. If the solver detects signs of stiffness—such as instability or rapidly shrinking step sizes—it switches to Backward Differentiation Formula (BDF) methods, which are implicit and more stable for stiff systems but computationally more expensive because they require solving nonlinear equations at each step. Throughout the integration, LSODA continuously adjusts the step size and method order to satisfy user-specified error tolerances, typically expressed in terms of relative and absolute error bounds. By combining automatic stiffness detection with adaptive step size and order control, LSODA provides a robust solver that performs well across a wide range of ODE problems without requiring detailed tuning from the user. It is likewise simple to implement with solve_ivp.

```python
from scipy.integrate import solve_ivp

def LSODA(X0, Y0, tmin, tmax, nts, du_dt):
    t_span = (tmin,tmax)
    t = np.linspace(tmin,tmax,nts,endpoint = False)

    #Express ICs as arrays
    X0 = np.asarray(X0)
    Y0 = np.asarray(Y0)

    #Again, repackaging the solutions
    solutions = [t]
    for i in range(len(X0)):
        sol = solve_ivp(du_dt, t_span, [X0[i], Y0[i]], t_eval = t, method = 'LSODA')
        solutions.append([sol.y[0], sol.y[1]])

    return solutions
```
The last two algorithms belong to a class of methods known as symplectic integrators. When possible, they conserve the volume for a continuous patch of initial condtions in phase space as it evolves in time. For physics problems, this translates to conseravation of energy. This makes them useful for long-term simulations of systems that conserve energy. 

The first one to be considered is the velocity Verlet method. It is a second order symplectic method that works by updating position and velocity in a way that incorporates how acceleration changes over time instead of treating it as constant over a full step. Intuitively, you first use the current velocity and acceleration to “predict” where the particle will move over a small time step, giving a very accurate new position. Then, because forces may depend on position, you recompute the acceleration at this new location. Finally, you update the velocity using the average of the old and new accelerations, which captures how the force changed during the step. the solve_ivp function does not have this method as an option, so it has to be written explicitly.

```python
def VelVerlet(X0, Y0, tmin, tmax, nts, du_dt):
    t = np.linspace(tmin,tmax,nts,endpoint = False)

    #Doing this allows me to vectorize the solver, so X0 and Y0 can be arrays of initial conditions, and allows coupling.
    X0 = np.asarray(X0)
    Y0 = np.asarray(Y0)

    #Make X and Y the same  structure as the 
    X = np.zeros((len(t), ) + X0.shape)
    Y = np.zeros((len(t), ) + Y0.shape)

    dt = t[1]-t[0]

    X[0] = X0
    Y[0] = Y0

    for it in range(0,nts-1):
        X[it+1] = X[it] + dt*Y[it] + 0.5*(dt**2)*(du_dt(t[it],[X[it], Y[it]])[1])

        Y_predict = Y[it] + dt*(du_dt(t[it],[X[it], Y[it]])[1])

        Y[it+1] = Y[it] + 0.5*dt*(du_dt(t[it],[X[it], Y[it]])[1] +  du_dt(t[it+1], [X[it+1],Y_predict])[1])

    # Now, we need to repackage these results so that they are implemented the same as the scipy methods
    if len(X0) == 1: #special case for scalar inputs
        X = [X[i][0] for i in range(len(X))]
        Y = [Y[i][0] for i in range(len(Y))]
        return [t, [np.asarray(X), np.asarray(Y)]]
    else:
        solutions = [t]
        for k in range(len(X0)):
            U = [X[i][k] for i in range(len(X))]
            V = [Y[i][k] for i in range(len(Y))]
            solutions.append([np.asarray(U), np.asarray(V)])

        return solutions
```
The last method that will be used is the 4th Order Yoshida Integrator. It was devloped by plasma physicists before making its way into the other fields. The Yoshida method builds a very accurate time step by carefully combining several smaller, symmetrically arranged steps of a simpler symplectic integrator (like velocity Verlet). Instead of taking one step with a fixed update rule, it takes a sequence of substeps with specially chosen positive and negative time coefficients that cause lower-order errors to cancel out. The result is a method that remains symplectic and achieves much higher accuracy per step than other methods. Likewise, this will have to be written explicitly.

```python
def Yoshida(X0, Y0, tmin, tmax, nts, du_dt):
    t = np.linspace(tmin,tmax,nts,endpoint = False)

    X0 = np.asarray(X0)
    Y0 = np.asarray(Y0)

    X = np.zeros((len(t), ) + X0.shape)
    Y = np.zeros((len(t), ) + Y0.shape)

    dt = t[1]-t[0]

    X[0] = X0
    Y[0] = Y0

    #Yoshida coefficients
    w1 = 1/(2 - np.cbrt(2))
    w2 =  -np.cbrt(2)/(2 - np.cbrt(2))

    c1 = w1/2
    c2 = (w1 + w2)/2
    c3 = c2
    c4 = c1

    d1 = w1
    d2 = w2
    d3 = w1
    
    #Looping to do the method
    for it in range(0,nts-1):
        Y1 = Y[it] + c1*dt*(du_dt(t[it], [X[it], Y[it]])[1])
        X1 = X[it] + d1*dt*(du_dt(t[it], [X[it], Y1])[0])
        t1 = t[it] + d1*dt

        Y2 = Y1 + c2*dt*(du_dt(t1, [X1, Y1])[1])
        X2 = X1 + d2*dt*(du_dt(t1, [X1, Y2])[0])
        t2 = t1 + d2*dt

        Y3 = Y2 + c3*dt*(du_dt(t2, [X2, Y2])[1])
        X3 = X2 + d3*dt*(du_dt(t2, [X2, Y3])[0])
        t3 = t2 + d3*dt

        Y[it+1] = Y3 + c4*dt*(du_dt(t3, [X3, Y3])[1])
        X[it+1] = X3

    # Now, we need to repackage these results so that they are structured the same as the scipy methods
    if len(X0) == 1:
        X = [X[i][0] for i in range(len(X))]
        Y = [Y[i][0] for i in range(len(Y))]
        return [t, [np.asarray(X), np.asarray(Y)]]
    else:
        solutions = [t]
        for k in range(len(X0)):
            U = [X[i][k] for i in range(len(X))]
            V = [Y[i][k] for i in range(len(Y))]
            solutions.append([np.asarray(U), np.asarray(V)])

        return solutions
```
One last thing to note is the limits of these function. The way the RK4(5) and LSODA funtions are written, it is not possible to have a system where particles are couple together. This is because the functions merely loop over intial condtions and solve that particular problem. Since the the symplectic methods had to be written explicitly, they have much more freedom in what they can do. If you write your system correctly and express the derivative properly, these methods can treat the lists of initial conditions as interacting particles. Of course, if you instead want to solve the same problem for different sets of inital conditions, that is still possible too. You just don't couple your derivatives of one trajectory to the others. This makes them highly flexible! They take full advantage of the vectorization of numpy arrays. The symplectic methods can solve something like the SHO for multiple ICS up to the n-body problem for as many bodies as you input without changing the function at all. We will see all of these uses in the next sections.

## Simple Harmonic Oscillator

The first system that we will investigate the simple harmonic oscillator with and without damping. To simplify the equations, the equillbrium position is taken to be the origin so that erroneous terms accounting for the length of the spring do not have to be accounted for. Newton's Second Law gives the following equations of motion.

```math
\ddot{x} = \frac{c}{m} \dot{x} + \frac{k}{m} x
```
Here, k is the spring constant, m is the mass, and c is the damping strength. For the integrators, we need to express this system as a coupled system of first order equations. 

```math
\dot{x} = y
```
```math
\dot{y} = -\frac{b}{m} y - \frac{k}{m} x
```
For the way my code is configured, we need to express the derivatives of the system in a specfic way.
```python

k = 1 # N/m        (spring constant)
m = 1 # kg         (mass)
c = 0.5 # N/(m/s)  (damping strength)

def SHO_damped(t, u): #derivative, u = [X, Y]
    x,y = u

    #Put the derivative for x and y in here
    #x is the position. y is the velocity
    dx_dt = y
    dy_dt = -(k/m)*x  - (c/m)*y
    return [dx_dt, dy_dt]
```
Let's first look at how the this system evolves over time with no damping for the different integrators. 

<div align="center">
  <img src="SHO_Phase.png" alt="Undamped Trajectories" width="600">
  <p><em>Figure 1:</em> Plots of the phasae space trajecotries for the undamped SHO for three different integrators.</p>
</div>

We can see that this traces out an ellipse in phase space, which is what we expect for the analytic solution for system. However, these methods are not outputting the exact same trajectory in phase space. You can resolve the difference by zooming in the plot range, but this is a good oppurtunitty to see how the total mechanical enerergy of this system evolves with the different methods. Physically, we know that it should be conserved because spring forces are conservative, but not all numeriacl methods conserve energy.

<div align="center">
  <img src="SHO_Energy.png" alt="Undamped Energies" width="600">
  <p><em>Figure 2:</em> Plots of the energies for the undamped SHO for three different integrators.</p>
</div>

This shows a large distinction among these three methods. We see that the symplectic method oscillates near the true energy, while the RK4(5) and LSODA methods accumulate error over time that causes them to drift from the true solution. Since RK4(5) and LSODA unavoidably accumulate drift in the energy over time, Veleocity Verlet might be considered better for long term solutions to this problem despite being a lower order method. Of course, this isn't a complicated problem, but, for somehing more complex, energy drift is signifigant. Now, we can introduce damping to this system and see how that affects the numeriacl solutions. 

<div align="center">
  <img src="SHO_Damped_Phase.png" alt="Damped Trajectories" width="600">
  <p><em>Figure 3:</em> Plots of the phase space trajectories for the damped SHO for three different integrators.</p>
</div>

As expected, the system is now losing energy to viscous forces. This causes it to spiral down to the origin, which for this system is the state with no energy. We see the same behavior for the energy of this system.

<div align="center">
  <img src="SHO_Damped_Energy.png" alt="Damped Energies" width="600">
  <p><em>Figure 4:</em> Plots of the energy over time for the damped SHO for three different integrators.</p>
</div>

In this case, Velocity Verlet just behaves as a second order integrator. Therefore, it has lost the advantage it had over RK4(5) and LSODA.

## n - Body Problem and the Virial Theorem
Too simulate the n-body problem, we need to develop some more complicated theory than the SHO. We can derive approximate equations for the gravitational field produced by a point of mass M by using the weak-field limit of the Einstein field equations. When doing so, we can obtain Newton's universal law of gravitation. 

```math
\mathbf{g}(\mathbf{x}) = -G M \frac{\mathbf{x}}{x^3}
```
Where x is the position vector from the mass point to the field point, and G is the gravitational constant. Then, for a collection of n point masses, the acceleration on each point is the sum of the individual accelerations. 

```math
\mathbf{a}_j = - \sum_{\substack{i=1 \\ i \neq j}}^{n} G m_i \frac{\mathbf{x}_j - \mathbf{x}_i}{\lvert \mathbf{x}_j - \mathbf{x}_i \rvert^3}
```
Now, the positon vectors are taken to be relative to the origin of the system. We want to express this a coupled first order system, so let y denote velocities.
```math
\dot{\mathbf{x}}_j = \mathbf{y}_j,
\qquad
\dot{\mathbf{y}}_j = - \sum_{\substack{i=1 \\ i \neq j}}^{n} G m_i \frac{\mathbf{x}_j - \mathbf{x}_i}{\lvert \mathbf{x}_j - \mathbf{x}_i \rvert^3}
```
Now, we've essentially written this as a first order system. We can repackage these results to resemble oher system. Let X be the vector whose components are the positions x_j, then let Y be the vector whose components are y_j, lastly let A be the vector whose components are a_j. Then, our system can be written in the highly absctraced form

```math
\frac{d}{dt}
\begin{pmatrix}
\mathbf{X} \\
\mathbf{Y}
\end{pmatrix}
=
\begin{pmatrix}
\mathbf{Y} \\
\mathbf{A}
\end{pmatrix}
```
This may not seem useful, but it means that we just need to calculate A as a list of 2D vectors at each timestep to simulate this problem. We can do that with the symplectic integrators because of how general they are due to numpy arrays being vectorized. To compute the derivatives, we use the following functions. This works well anyways since symplectic methods are especially useful for n-body simulation.

```python
def GravAcc(P1, P2, M2): #acceleration on P1 from P2
    if np.array_equal(P1, P2): #Impportant edge case
        return np.array([0,0])
    else:
        return G*M2*(P2 - P1)/((np.linalg.norm(P2 - P1)**3)) #Gravitational acceleration formula
    
def TotalAcceleration(P, X, M): # Calculates acceleratio on P from all points in X
    A = np.array([0, 0])
    for i in range(len(X)):
        A = A + GravAcc(P, X[i], M[i])
    return A

def Derivative(t, U): # U is a 2 element list. The first element is the array of position. The second is velocities
    X, Y = U

    dX_dt = Y
    dY_dt = np.asarray([TotalAcceleration(X[i], X, M) for i in range(len(X))])
    return [dX_dt, dY_dt]
```

Then, we can test this with five points masses. We are able to write our conditions for the system as follows.

```python
G = 1 #Gravitational Constant 
M =  [1.0, .50, 0.50, .450, 0.450]            # Masses
X0 = [[0,0], [2,0], [-2,0], [0,2], [0,-2]]    # Initial Positions
V0 = [[0,0], [0,1], [0,-1], [-1,0], [1,0]]    # Initial Velocities
```
Then, the simulation results are plotted.

<div align="center">
  <img src="5_Body_Trajectories.png" alt="5-Body simulation" width="600">
  <p><em>Figure 5:</em> Trajectories of 5 gravitationally attracting masses.</p>
</div>

We have enough particles in this simulation to test the virial theorem. For a gravtiationally bound system, it related the average kinetic and potential energy.

```math
\langle T \rangle = -\frac{1}{2} \langle U \rangle
```

Let's first see how the kinetic and potential energy evolves over time.

<div align="center">
  <img src="5_Body_Energies.png" alt="5-Body Energies" width="600">
  <p><em>Figure 6:</em> total potential and kinetic energy of for 5 gravitationally attracting masses.</p>
</div>

The energy is clearly being conserved, and when we calculate the averages of kinetic and potential, we find the following.

```math
\langle T \rangle = 0.5790...
```
```math
\langle U \rangle =-1.0103...
```
When we compute their ratio, we find
```math
\frac{\langle T \rangle}{\langle U \rangle} = -0.5249...
```
Which is pretty close to the -1/2 we would expect, so we can be reasonably sure that this is confirming the virial theorem.

## Phase Space Volume and Symplectic Methods

Earlier, symplectic methods were described as conserving volume in phase space. This is something we can analyze with our undamped simple harmonic oscillator. Since this system follows the trajectory of one particle, we have a two dimensional phase space: Velocity vs Position. So, our volume becomes an area. Then, the algorithm is relativly simple. We want to see how the volume of a small patch of area in phase space evolves over time. To achieve this, we generate a list of random initial conditions that are close. This will approximate our patch of area. Then, we will use the concave_hull and shapely.geometry packages to generate a concave hull around our points as they evolve in time and calculate the area. Here is some sample code for the Velocity Verlet method.

```python
from concave_hull import concave_hull
from shapely.geometry import Polygon
from P3_Methods import VelVerlet, Yoshida, LSODA

tmin = 0 #s start time
tmax = 50 #s end time
nts = 350 #number of points between tmin and tmax

n = 250
X0 = np.random.uniform(0.95, 1, size = n)
Y0 = np.random.uniform(0.95, 1, size = n)

solutions = VelVerlet(X0, Y0, tmin, tmax, nts, SHO)
t = solutions[0]

PhaseVolumes = np.zeros(len(t))
for j in range(len(t)):
    P = [[solutions[k][0][j], solutions[k][1][j]] for k in range(1,n)]
    P = concave_hull(P, concavity=2.0)
    poly = Polygon(P)
    PhaseVolumes[j] = poly.area

plt.plot(t, PhaseVolumes, label = "Velocity Verlet")
```
We can see how these different methods affect phase space. 

<div align="center">
  <img src="SHO_Volume.png" alt="Undamped Phase Volumes" width="600">
  <p><em>Figure 7:</em> Plots of the phasae space volums over time for the undamped SHO for three different integrators.</p>
</div>

The two symplectic methods seem to be conserving phase space volume while the others are accumulating more error as the simulation continues! Now, let's see what happens when damping is introduced.

<div align="center">
  <img src="SHO_Damped_Volume.png" alt="Damped Phase Volumes" width="600">
  <p><em>Figure 8:</em> Plots of the phasae space volums over time for the damped SHO for three different integrators.</p>
</div>

Now, all methods have the volume decaying to zero. This is not surprising of course. This system no longer conserves energy, and all initial conditions should tend towards the zero energy state at the origin. So, since all points get closer over time, the phase space area will decrease.

## Conclusion
(Idea to inclulde in this section: Which integrator is the best? There is not one that is objectively better than the others. The one that you choose depends on the the problem you are trying to solve and the context of your field.)

## Attribution
The book 

## Timekeeping
I spent about 20 hours on the code and about 10 hours on the write up.

## Languages, Libraries, Lessons Learned
Everything was written in python. I made an effort to write all of my code to be highly reusable. This meant that the same function I use to simulate the SHO could be used to simulate the n-body problem. 
