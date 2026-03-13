# Project 3
## Introduction
Few areas of math enjoy such a privelleged position in physics as differential equations do. However, for all of their utility in modeling physical system, the average differential equation(s) that we encounter outisde of an introductory course lack any kind of analytic solution. This necessitates the development of numerical methods that can efficiently and accuratley approximate solutions to initial value problems. In this porject, we will examine different algorithms for solving ODEs and see how they solve the simple harmonic oscillator along with some other special cases.

## Algorithms and Theory
We will first examine the different numerical methods that will be used. The first one is RK4(5). It belongs to a family of solvers known as Runge-Kutta (RK) methods. The explicit derivation is not of interest here, but a brief explanation of how this family of solvers works is useful. They are what's known as prediction corrector methods. Unlike bad solvers such as Euler's method, "predictor-corrector methods improve the approximation accuracy by querying the 𝐹 function several times at different locations (predictions), and then using a weighted average of the results (corrections) to update the state." (INJECT THE REST OF AN OUTLINE FOR THE DERIVATION HERE). One last note is that RK4(5) is unique in that it is an adaptive method. It achieves by chagning the step size through comparing a 4th order step and a 5th order step. A python implementation is simple using premade librairies. 

```python
from scipy.integrate import solve_ivp

def RK45(X0,Y0,tmin,tmax,nts,du_dt):
    t_span = (tmin,tmax)
    t = np.linspace(tmin,tmax,nts,endpoint = False)
    U0 = [X0,Y0]

    return solve_ivp(du_dt, t_span, U0,t_eval = t, method = 'RK45')
```

The next method used is known as LSODA. It is an adaptive numerical solver for systems of ordinary differential equations. It was developed as part of the ODEPACK library and is designed to efficiently handle both stiff and non-stiff problems without requiring the user to decide which type of solver to use. It automatically switches between two classes of multistep methods depending on the behavior of the system. When the problem appears non-stiff, LSODA uses variable-order Adams predictor–corrector methods, which are explicit multistep schemes that are efficient for smooth solutions. If the solver detects signs of stiffness—such as instability or rapidly shrinking step sizes—it switches to Backward Differentiation Formula (BDF) methods, which are implicit and more stable for stiff systems but computationally more expensive because they require solving nonlinear equations at each step. Throughout the integration, LSODA continuously adjusts the step size and method order to satisfy user-specified error tolerances, typically expressed in terms of relative and absolute error bounds. By combining automatic stiffness detection with adaptive step size and order control, LSODA provides a robust solver that performs well across a wide range of ODE problems without requiring detailed tuning from the user. It is likewise simple to implement with premade libraries.

```python
from scipy.integrate import solve_ivp

def LSODA(X0,Y0,tmin,tmax,nts,du_dt):
    t_span = (tmin,tmax)
    t = np.linspace(tmin,tmax,nts,endpoint = False)
    U0 = [X0,Y0]

    return solve_ivp(du_dt, t_span, U0,t_eval = t, method = 'LSODA')
```
The last algorithm to discuss is called Velocity Verlet. (Put derivation and explanation for velocity verlet here). Here is the implementation for python.
```python
def VelVerlet(X0, Y0, tmin, tmax, nts, du_dt):
    t = np.linspace(tmin,tmax,nts,endpoint = False)
    X = np.zeros(nts)
    Y = np.zeros(nts)
    U = np.empty(nts, dtype=object)

    dt = t[1]-t[0]

    X[0] = X0
    Y[0] = Y0
    U[0] = [X0,Y0]

    for it in range(0,nts-1):
        X[it+1] = X[it] + dt*Y[it] + 0.5*(dt**2)*(du_dt(t[it],U[it])[1])

        Y_predict = Y[it] + dt*(du_dt(t[it],U[it])[1])

        Y[it+1] = Y[it] + 0.5*dt*(du_dt(t[it],U[it])[1] +  du_dt(t[it+1], [X[it+1],Y_predict])[1])

        U[it+1] =[X[it+1],Y[it+1]]


    return t, X, Y
```

## Applications and Discussion

The first system that we will investigate the simple harmonic oscillator with and without damping. To simplify the equations, the equillbrium position is taken to be the origin so that erroneous terms for the length of the spring do not have to be accounted for. (Talk about the functions used to compute the derivatives for the system, put the code here, then show the plots for the unddamped and damped case. Proceed to talk about how the verlet method conserves energy and does not drift, and how to other methods maintain higher precision when the energy is not conserved. Don;t forget to discuss the order of error and how that impacts the long-term accuracy of the method.)

## Conclusion

## Attribution

## Timekeeping

## Languages, Libraries, Lessons Learned
