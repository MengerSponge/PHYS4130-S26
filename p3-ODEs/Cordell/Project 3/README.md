# Introduction
Numerical methods can be implemented to approximate solutions for ordinary differential equations (ODEs) that would otherwise be a computationally labor-some task. This package was developed to demonstrate the capabilities of several approximation techniques which include Euler's method, the Runge-Kutta technique, Verlet integration, and Scipy's ODEINT. We provide phase-space and Energy vs Time plots to illustrate their approximation capabilities and ability to conserve energy. By the same token, analyze the absolute error of each method.        


 We found that the Verlet method reaches a 5% relative error within 64 time steps, whereas RK2, RK4, and ODEINT reach the target error within 128 time steps, and Euler's method requires about 2048 time steps. The pros and cons of each method depending on the situation are further discussed in the report.
## Background Theory

Euler's method is typically what one would start out with when exploring ODE approximation methods. It approximates the solution of an ODE at a point, B, by starting from an initial value, point A, and taking the next point B on the tangent line to the solution at point A. Then, it repeats the process for the subsequent points as depicted in Fig. 1. The number of these sub-intervals taken over time is typically referred to as the number of time steps (nts). And, the error of Euler's method reduces proportionally as the number of time steps increases, so we call it a first order method or a first order approximation. 

<p align="center">
  <img src="./Euler_method.png" alt="Extension plot" width="300">
</p>

<p align="center">
  Figure 2: Illustration of Euler's method. The red line is the numerical approximation, and the blue line is the analytic solution. Reproduced from [1].
</p>



The Runge-Kutta method is similar to Euler's method in the sense that it is an iterative technique; however, it incorporates additional points within the interval which generally makes it much more accurate then Euler's method. For example, RK2 uses the slope at A to approximate the solution at the midpoint of A and B, then uses the slope at the midpoint to approximate the solution at B. Furthermore, RK4 utilizes a weighted average of the slope at A, the slope at B, and the slope at the midpoint using two different estimates. As a result, it turns out that RK2 is a second order method whereas RK4 is a fourth order method.

Verlet integration can be implemented for second order ODE's of the form $\ddot{x}(t) = A(x(t))$ such as for a harmonic oscillator. The algorithm is derived from the Taylor expansion for $x(t+\Delta t)$ and $x(t-\Delta t)$ as follows,

$$
x(t+\Delta t)=x(t)+\dot{x}(t)\Delta t+\frac{1}{2}\ddot{x}(t)\Delta t^2+\frac{1}{6}x^{(3)}(t)\Delta t^3+\cdots
$$

$$
x(t-\Delta t)=x(t)-\dot{x}(t)\Delta t+\frac{1}{2}\ddot{x}(t)\Delta t^2-\frac{1}{6}x^{(3)}(t)\Delta t^3+\cdots
$$

where adding these yields

$$
x(t+\Delta t)=2x(t)-x(t-\Delta t)+\ddot{x}(t)\Delta t^2+\mathcal{O}(\Delta t^4)
$$

Scipy's ODEINT incorporates a bunch of stuff it seems...

# Procedure
To demonstrate the use of each method, we approximate x(t) and v(t) for a harmonic oscillator with and without a linear dampening term. We illustrate these solutions by plotting the phase-space diagrams for each method. In addition, we provide Energy vs. Time plots to show whether or not a specific method conserves energy. To compare between the different methods, our tool offers several means of error analysis where the error is computed from the analytic solution. The first is a straightforward output of the error at some user-selected time. This option is most reliable (as it works for RK4 and ODEINT), yet it offers the least insight. The next option is to output the number of time steps required to achieve a user-defined target error. However, this option is not compatible with RK4 and ODEINT as Scipy adaptively chooses the number of timesteps. Lastly, our tool offers Number of Time Steps vs Error loglog plots which also is not compatible with RK4 and ODEINT since Scipy adaptively chooses the number of timesteps. Initially, we attempted to work around Scipy's adaptive selection of nts by passing specific time points in an array, t_eval, "linspaced" into time points based on the number of time steps. However, Scipy just returns the solutions at those time points based on its own number of time steps. 

 
## Function Code

### Euler_Solver

```python
def SHO_solver_Euler(x0, v0, tmin, tmax, nts, SHO_deriv):
    x_array = np.zeros(nts)                           
    v_array = np.zeros(nts)                           
    t_array = np.linspace(tmin, tmax, nts, endpoint=False)
    dt = t_array[1] - t_array[0]                        
    x_array[0] = x0                                     
    v_array[0] = v0                                     
    
    for it in range(0, nts-1):
        x_array[it+1] = x_array[it] + dt * SHO_deriv([x_array[it], v_array[it]], t_array[it])[0]
        v_array[it+1] = v_array[it] + dt * SHO_deriv([x_array[it], v_array[it]], t_array[it])[1]
    
    return t_array, x_array, v_array
```
The recursive algorithm we use for Euler's method shown above is: 
$$
x_0 = x_0 \\
v_0 = v_0
$$

$$
x_{n+1} = x_n + \Delta t\,x_n   \\
v_{n+1} = v_n + \Delta t\,v_n
$$

### RK2_Solver
```python
def SHO_solver_RK2(x0, v0, tmin, tmax, nts, SHO_deriv):
    x_array = np.zeros(nts)                                                     # array to hold position
    v_array = np.zeros(nts)                                                     # array to hold velocity
    t_array = np.linspace(tmin, tmax, nts, endpoint=False)                      # array holds the time points 
    dt = t_array[1] - t_array[0]                                                # dt = time step length  
    x_array[0] = x0                                                             # Initial position
    v_array[0] = v0                                                             # Initial velocity
    for it in range(0, len(t_array)-1 ):                                        # loop over time steps
        t  = t_array[it]                                                        
        x_h = x_array[it] + (dt/2 * SHO_deriv([x_array[it], v_array[it]], t)[0])                # sub-step 1 for RK2
        v_h = v_array[it] + (dt/2 * SHO_deriv([x_array[it], v_array[it]], t)[1])                # sub-step 1 for RK2
        x_array[it+1] = x_array[it] + (dt * SHO_deriv([x_h, v_h], t + dt/2)[0])         # sub-step 2 for RK2
        v_array[it+1] = v_array[it] + (dt * SHO_deriv([x_h, v_h], t + dt/2)[1])         # sub-step 2 for RK2
    return t_array, x_array, v_array
```
The recursively defined algorithm for RK2_solver is:
   $$
   x_{n+\frac{1}{2}} = x_n + \frac{\Delta t}{2}\,v_n
   $$

   $$
   v_{n+\frac{1}{2}} = v_n + \frac{\Delta t}{2}\,A(x_n,v_n)
   $$

   
   $$
   x_{n+1} = x_n + \Delta t\,v_{n+\frac{1}{2}}
   $$

   $$
   v_{n+1} = v_n + \Delta t\,A\!\left(x_{n+\frac{1}{2}}, v_{n+\frac{1}{2}}\right)
   $$
   where A is the acceleration.
### Verlet_Solver

```python
def verlet_solver(x0, v0, tmin, tmax, nts, deriv):
    x_array = np.zeros(nts)                                                     # array to hold position
    v_array = np.zeros(nts)                                                     # array to hold velocity                                               
    t_array = np.linspace(tmin, tmax, nts, endpoint=False)                      # array holds the time points 
    dt = t_array[1] - t_array[0]                                                # dt = time step length  
    x_array[0] = x0                                                             # Initial position
    v_array[0] = v0                                                             # Initial velocity
    for it in range(0, len(t_array)-1 ):                                        # loop over time steps
        # Algorithm for Verlet method 
        x1 = x_array[it] + v_array[it]*dt + 0.5 * deriv(x_array[it], v_array[it]) * dt**2
        v1 = v_array[it] + 0.5 * (deriv(x_array[it], v_array[it]) + deriv(x1, v_array[it])) * dt
        x_array[it+1] = x1
        v_array[it+1] = v1


    return t_array, x_array, v_array
```
The verlet algorithm above is recursively defined as:


   $$
   x_0 = x_0 \\
    v_0 = v_0
   $$

   $$
   x_{n+1} = x_n + v_n\,\Delta t + \frac{1}{2}A(x_n,v_n)\,\Delta t^2
   $$

   $$
   v_{n+1} = v_n + \frac{1}{2}\left[A(x_n,v_n) + A(x_{n+1},v_n)\right]\Delta t
   $$
where A is acceleration.


## Instructions
To run the tool use this command from the terminal that includes P3_code.py and Functions.py.

```bash
conda activate [environment]
python P3_Code.py
```



# Analysis



## Phase Space 



## Energy vs. Time 



## Error 



# Conclusions
This is where I want to talk about in what cases might one want to use which method like in terms of computational cost and overall efficiency.

# Extensions

# Questions

## Timekeeping

Week before spring break: 
1 hour: Tuesday in class

1 hour: Tuesday after class or Wednesday (I forgot which day)

1? hour: Thursday

1 hour: Friday

After spring break: 
3 hours: Tuesday 3/24
