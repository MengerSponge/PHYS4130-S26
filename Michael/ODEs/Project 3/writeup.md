# Project 3
## Introduction
Few areas of math enjoy such a privelleged position in physics as differential equations do. However, for all of their utility in modeling physical system, the average differential equation(s) that we encounter outisde of an introductory course lack any kind of analytic solution. This necessitates the development of numerical methods that can efficiently and accuratley approximate solutions to initial value problems. In this porject, we will examine different algorithms for solving ODEs and see how they solve the simple harmonic oscillator along with some other special cases.

## Algorithms and Theory
We will first examine the different numerical methods that will be used. The first one is RK4(5). It belongs to a family of solvers known as Runge-Kutta (RK) methods. The explicit derivation is not of interest here, but a brief explanation of how this family of solvers works is useful. They are what's known as prediction corrector methods. Unlike bad solvers such as Euler's method, "predictor-corrector methods improve the approximation accuracy by querying the 𝐹 function several times at different locations (predictions), and then using a weighted average of the results (corrections) to update the state." (INJECT THE REST OF AN OUTLINE FOR THE DERIVATION HERE). One last note is that RK4(5) is unique in that it is an adaptive method. It achieves by chagning the step size through comparing a 4th order step and a 5th order step. A python implementation is simple using premade librairies. 

```python
from scipy.integrate import solve_ivp

def RK45(X0, Y0, tmin, tmax, nts, du_dt):
    t_span = (tmin,tmax)
    t = np.linspace(tmin,tmax,nts,endpoint = False)
    
    X0 = np.asarray(X0)
    Y0 = np.asarray(Y0)

    solutions = [t]
    for i in range(len(X0)):
        sol = solve_ivp(du_dt, t_span, [X0[i], Y0[i]], t_eval = t, method = 'RK45')
        solutions.append([sol.y[0], sol.y[1]])

    return solutions
```

The next method used is known as LSODA. It is an adaptive numerical solver for systems of ordinary differential equations. It was developed as part of the ODEPACK library and is designed to efficiently handle both stiff and non-stiff problems without requiring the user to decide which type of solver to use. It automatically switches between two classes of multistep methods depending on the behavior of the system. When the problem appears non-stiff, LSODA uses variable-order Adams predictor–corrector methods, which are explicit multistep schemes that are efficient for smooth solutions. If the solver detects signs of stiffness—such as instability or rapidly shrinking step sizes—it switches to Backward Differentiation Formula (BDF) methods, which are implicit and more stable for stiff systems but computationally more expensive because they require solving nonlinear equations at each step. Throughout the integration, LSODA continuously adjusts the step size and method order to satisfy user-specified error tolerances, typically expressed in terms of relative and absolute error bounds. By combining automatic stiffness detection with adaptive step size and order control, LSODA provides a robust solver that performs well across a wide range of ODE problems without requiring detailed tuning from the user. It is likewise simple to implement with premade libraries.

```python
from scipy.integrate import solve_ivp

def LSODA(X0, Y0, tmin, tmax, nts, du_dt):
    t_span = (tmin,tmax)
    t = np.linspace(tmin,tmax,nts,endpoint = False)
    
    X0 = np.asarray(X0)
    Y0 = np.asarray(Y0)

    solutions = [t]
    for i in range(len(X0)):
        sol = solve_ivp(du_dt, t_span, [X0[i], Y0[i]], t_eval = t, method = 'LSODA')
        solutions.append([sol.y[0], sol.y[1]])

    return solutions
```
The last two algorithms to discuss fall under a class of integrators known as symplectic methods. They conserve the volume for a continuous patch of initial condtions in phase space that are evolved in time. For physics problems, this translates to conseravation of energy. This makes them useful for simulating systems that conserve energy over long time periods. The first one to be considered is the velocity Verlet method. It is a second order symplectic method that works by updating position and velocity in a way that naturally incorporates how acceleration changes over time instead of treating it as constant over a full step. Intuitively, you first use the current velocity and acceleration to “predict” where the particle will move over a small time step, giving a very accurate new position. Then, because forces (and thus acceleration) may depend on position, you recompute the acceleration at this new location. Finally, you update the velocity using the average of the old and new accelerations, which captures how the force changed during the step. This “average acceleration” idea is what makes velocity Verlet both more accurate and better at conserving energy over long simulations compared to simpler methods like Euler.

```python
def VelVerlet(X0, Y0, tmin, tmax, nts, du_dt):
    t = np.linspace(tmin,tmax,nts,endpoint = False)

    #Doing this allows me to vectorize the solver, so X0 and Y0 can be arrays of initial conditions.
    X0 = np.asarray(X0)
    Y0 = np.asarray(Y0)

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
The last method that will be used is the 4th Order Yoshida Integrator. It was historically devloped by plasma physicists before making its way into the other fields. The Yoshida method builds a very accurate time step by carefully combining several smaller, symmetrically arranged steps of a simpler symplectic integrator (like velocity Verlet). Instead of taking one step with a fixed update rule, it takes a sequence of substeps with specially chosen positive and negative time coefficients that cause lower-order errors to cancel out. Intuitively, you can think of it as moving the system forward, slightly backward, and forward again in a balanced way so that the mistakes made in each part offset one another. The result is a method that remains symplectic (so it preserves the geometric structure of the system, like energy behavior over long times) but achieves much higher accuracy per step than basic methods.

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
## Simple Harmonic Oscillator

The first system that we will investigate the simple harmonic oscillator with and without damping. To simplify the equations, the equillbrium position is taken to be the origin so that erroneous terms accounting for the length of the spring do not have to be accounted for. 

(LaTex to depict the ODE for the system, then show the coupled 1st order system)

For the way my code is configured, we need to express the derivatives of the system in a specfic way.
```python

k = 1 # N/m        (spring constant)
m = 1 # kg         (mass)
c = 0.5 # N/(m/s)  (damping term strength)

def SHO_damped(t, u): 
    x,y = u

    #Put the derivative for x and y in here
    #x is the position. y is the velocity
    dx_dt = y
    dy_dt = -(k/m)*x  - (c/m)*y
    return [dx_dt, dy_dt]
```
Let's first look at how the this system evolves over time with no damping for the different integrators. 

(Show the plot for the undamped oscillator here. Pick a time step good enough to show the different order errors active here)

We can see that this traces out an ellipse in phase space, which is what we expect for the analytic solution for system. However, these methods are not outputting the exact same trajectory in phase space. You can resolve the difference by zooming in the plot range, but this is a good oppurtunitty to see how the total mechanical enerergy of this system evolves with the different methods. Physically, we know that it should be conserved because spring forces are conservative, but that is not exactly the case for these methods.

(show plot of energy vs time to compare the methods)

This shows a large distinction among the methods! We see that the symplectic method oscillates near the true energy, whle the RK4(5) and LSODA methods accumulate error over time that causes them to drift from the true solution. For that reason, the symplectic Veleocity Verlet might be considered better for this problem despite being a lower order method. Now, we can introduce damping to this system and see how that affects the numeriacl solutions. 

(Show the phase space plots for the damped system )

As expected, the system is now losing energy to viscous forces. This causes it to spiral down to the origin, which for this system is the state with no energy. We see the same kind of behavior for the energy of this system.

(Show the energy vs time plots)

Now, to see how the different methods behave, we should consider the relative error of the methods with the analytic solution. 



then show the plots for the unddamped and damped case. Proceed to talk about how the verlet method conserves energy and does not drift, and how to other methods maintain higher precision when the energy is not conserved. Don;t forget to discuss the order of error and how that impacts the long-term accuracy of the method.)

## n - Body Problem and the Virial Theorem

## Phase Space Volume Conservation

## Conclusion
(Idea to inclulde in this section: Which integrator is the best? There is not one that is objectively better than the others. The one that you choose depends on the the problem you are trying to solve and the context of your field.)

## Attribution

## Timekeeping

## Languages, Libraries, Lessons Learned
