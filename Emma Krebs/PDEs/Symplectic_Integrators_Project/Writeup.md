---
meta:
    author: Emma Krebs
    topic: Symplectic Integrations
    course: TN Tech PHYS4130
    term: Spring 2026
---

# Symplectic Integrations

## Types of Integrations

Ordinary Differential Equations (ODEs) are differential equations that depend on one independent variable. ODEs can relate how systems can change over time or space, making them an extremely pervasive tool of mathematics for understanding how something develops. These functions can be solved numerically through integrators given an initial state. Thus, it is important to understand what these integrators are doing to our system and how they might accumulate error. An accurate ODE solver will be one that conserves the energy of the system, but what does that mean in context of our problem? We will start by considering harmonic oscillation, which also happens to be the main system we will investigate throughout this report. The function is given as:

$$m\text{ }\dfrac{d^2y}{d^2t} = -ky$$

where its associated equation of energy is

$$E = \dfrac{1}{2}\text{ } m \text{ }(\text{ } \dfrac{dy}{dt} \text{ })^2 + \dfrac{1}{2}ky^2$$


### What is Harmonic Oscillation and What Does it Look Like?

The main system we studied with these methods of integration was harmonic oscillation, which is defined to be a system that experiences a restoring force proportional to its displacement towards an equilibirum point. Some classical examples of this are mass springs and pendulums, where the restoring forces are from the spring and gravity, respectively. However, the harmonic oscillator is even more important to physics than just describing these simple examples because any mass subject to a force in stable equilibrium acts as a harmonic oscillator in small vibrations. Therefore, there are many physical phenomenon we can apply this system to, such as: quantum's potential wells, orbital perturbations, stellar oscillations, and more. 

Harmonic oscillators are characterized by two main categories, damped and undamped oscillators. These change what their phase plot diagrams look like, which are their space vs momentum graphs. For undamped systems, they are circular or elliptical shaped paths depending on the angular frequency of the system. An angular frequency of one results in a circular shape. An angular frequency less than one stretches the circle in the positon axis' direction and above one stretches it into the momentum axis' direction. If an oscillator is damped, this means that there is a term that decreases the amplitude over time. In response, the phase space diagram breaks the circular/elliptical shape and spirals inward. An example of what these graphs look like are demonstrated below for the three integrators we will be investigating:

<img src="Images_ODE/PhaseSpaceSymp.png" width="325"> <img src="Images_ODE/PhaseSpaceOdeint.png" width="325"> <img src="Images_ODE/PhaseSpaceRK45.png" width="325">

*Fig 1. This figure demonstrates three phase space graphs. From left to right, we have the Symplectic, Odeint, and RK45 integration methods.*

Although these all maintain similar shapes, if you look closely some appear to have thicker lines. These are from an accumulation of errors due to the method of integration. Let us see what makes these integrations unique and what these errors look like!

### Symplectic Integration

- Explain what symplectic integration is and the function you used, what it does good at.
- Include image showing the process the function goes through with midpoints. Include equation of calculation it does.
- Explain why this makes sense.
- Talk about linear dampending (spiral!)

Symplectic is unique 

```python
    def Verlet_symplectic(x_0, p_0, tmax, w, damp, N):

        t_array = np.linspace(0, tmax, N)
        x_array = np.zeros(len(t_array))
        p_array = np.zeros(len(t_array))
        h = t_array[1] - t_array[0]
        x_array[0] = x_0
        p_array[0] = p_0
    
        # Get the next value
        a_0 = -w**2*x_0 - damp*p_0
        x_array[1] = x_0 + p_0*h + 0.5*a_0*h**2
    
        # Iterations, update momentum first then position
        for i in range(1, len(t_array) - 1):
    
            a = -w**2*x_array[i] - damp*((x_array[i] - x_array[i - 1])) / h # Acceleration
            x_array[i+1] = 2*x_array[i] - x_array[i - 1] + a*h**2
    
            p_array[i] = (x_array[i+1] - x_array[i - 1]) / (2*h)
        
        p_array[-1] = (x_array[-1] - x_array[-2]) / h
    
        return x_array, p_array, t_array
```


### RK45

RK45, also known as Runge–Kutta–Fehlberg method, is an adapative numerical technique for solving ODEs. It uses intermediate calculations to produce two estimates of different accuracy. By doing this, it can adjust the time step sizes for high accuracy and efficiency. These technqiues have a range of different orders you can use for a probelm. In this case, RK45 computes a 4th and 5th order error estimate while RK2 would be a lower order with only a 2nd order estimate. Although RK45 excels at accuracy, it does not maintain the geometrical properties of the phase space, so it loses total energy from small numerical errors over time.

The function definition in the code is shown as:

```python
    def RK45_solver(x_0, p_0, tmax, w, damp, N):

        y_0 = [x_0, p_0]
        t_array = np.linspace(0, tmax, N)
        
        # Force it to match other integration methods
        
        sol = solve_ivp(Harmonic_deriv, (0, tmax), y_0, method='RK45', t_eval=t_array, args=(w, damp)) 
        t_array = sol.t
        x_array = sol.y[0]
        p_array = sol.y[1]
    
        return x_array, p_array, t_array
```



- Include image showing the process the function goes through with midpoints. Include equation of calculation it does.
- Explain the snippet of code for Rk45 and how you set it up, what it does with its arguments, etc.

### Odeint

- Explain what the Odent function is and what it does better than the other 
- Explain the snippet of code for Odeint and how you set it up, what it does with its arguments, etc.
- Talk about linear dampending

Calls the odeint function from scipy.integrate. odeint solves a system of ordinary 

```python

    def Odeint_solver(x_0, p_0, tmax, w, damp, N):

        y_0 = [x_0, p_0]
    
        t = np.linspace(0, tmax, N)
        L = odeint(Harmonic_deriv, y_0, t, args=(w, damp), tfirst=True, rtol=1e-4)
        
        return L[:,0], L[:,1], t
```

This function definition calls the Odeint function from scipy.integrate, and Odeint solves a system of ordinary of ordinary 
differential equations using lsoda from the FORTRAN library odepack.

### Final Comparison

Now that we understand our three functions a little more, let us turn our attention to the main.py function. We saw before that for these integration methods that they had a pretty consistent phase space diagram. However, this was for a pretty small time step. This means that error isn't as easily accumulated between steps since not as much changes. This doesn't mean that there is no error, and if we compare it between two N spacing values N = 100 and N = 3000, we see how much this error can impact our energy conservation:

<img src="Images_ODE/PhaseSpace3000.png" width="400"> <img src="Images_ODE/PhaseSpaceN100.png" width="600"> 

*Fig 2. Here we have two phase space plots for our three integration methods and analytic solution. The time step is changed such that the spacing given by N is 3000 (on the left), resulting in small steps, and 100 (on the right), resulting in larger steps.*

Thus, the error is now apparent. However, we need a more direct way to compare it. To do this, we will be looking at the errors of these methods through Total Energy vs Time, Relative Error vs Time, and Energy Drift vs Time for a variety of N spacings. Let us start with the former. We need to start by finding the total energy using a function definition from ODE_Methods.py:

```python

    def Total_energy(x_array, p_array, w):

        kinetic_energy = []
        potential_energy = []
        total_energy = []
        
        # Assuming m = 1, we can find the kinetic, potential, and total energy through the following methods
        for x, p in zip(x_array, p_array):
            KE_value = (1/2) * p**2
            PE_value = (1/2) * w**2 * x**2 # Since w is angular frequency and m = 1, so k = w**2
            TE_value = KE_value + PE_value
    
            kinetic_energy.append(KE_value)
            potential_energy.append(PE_value)
            total_energy.append(TE_value)
        
        return kinetic_energy, potential_energy, total_energy
```

This function returns the kinetic, potential, and total energy for a given array of x and p (along with their angular frequency) for a harmonic oscillator. It assumes the mass equals one so it can easily calculate these energies. It uses the general equations for kinetic and potentital energy such that for each x and p in these arrays:

$$KE_{Value} = \dfrac{1}{2} * p^2$$
$$PE_{Value} = \dfrac{1}{2} * w^2 * x^2$$
$$TE_{Value} = KE_{Value} + PE_{Value} = \dfrac{1}{2} * p^2 + \dfrac{1}{2} * w^2 * x^2$$

where these values are then appended to their respective arrays and returned. Thus, we find that:

<img src="Images_ODE/TotalEnergyVSTime.png" width="800"> 

*Fig 3. Total energy vs Time for three integration methods. This graph spans over the time of 150 seconds, but was zoomed in between 0 to 25 seconds to easily see the trends between the energies.*

As we see above, our baseline of the analytic solution is a flat line, which makes sense since there should be no error in the harmonic oscillator function. Then, we see that the Symplectic and RK45 methods oscillate around the total energy value of the analytic solution. The Symplectic is expected because it conserves energy; however, the RK45 is decent at minimizing the accumulation of error due to its adjustments of its steps when integrating. We see near the end though that error starts accumulating enough to decrease the oscillation average, making RK45 not conserve energy. Finally, Odeint does not oscillate and is fairly sporadic with its total energy, eventually decreasing like RK45 but at a much quicker pace. 

To compare these differences further, let us look at Relative Error vs Time. Here we calculate the relative error for these energies using the analytic value and the equation:

$$\text{Relative Error} = \dfrac{|  \text{ analytic  -  estimated }  |}{(\text{analytic})}$$

and we find

<img src="Images_ODE/RelativeErrorVSTime.png" width="800"> 

*Fig 4. Relative error vs time graph for the three integration methods. This graph also spans over 150 seconds, but is limited to 25 seconds.*  

Again, we see the relative error oscillate between 0 and a positive value indicating that, even though there is an error, it is manageable and still allows conservation of energy. Odeint is obviously once again collecting relative error at an alarming rate and, although this graph doesn't extend past 25 seconds, RK45 is also steadily collecting error while oscillating at that error value. 

Finally, we will compare the energy drift over time for different step values. This will be the biggest indicator for stability and conservation for a method. 

<img src="Images_ODE/EnergyDriftVSTimeSymp.png" width="325"> <img src="Images_ODE/EnergyDriftVSTimeOdeint.png" width="325"> <img src="Images_ODE/EnergyDriftVSTimeRk45.png" width="325">

*Fig 5. Energy Drift vs Time for the three methods. Starting from the left, we have Symplectic, Odeint, and RK45.*

We see once again that even at low N the Symplectic method oscillates around 0, and for higher N it significantly decreases its oscillation range. Interestingly enough, the Symplectic is actually much worse than Odeint and RK45 for N = 100 if we look at the y-axis. Its range has a rough energy drift maximum of 2 while the others are on the order of $$10^{-3}$$ and $$10^{-4}$$. Although we know that Symplectic is overall the better choice for stability of harmonic oscillators, we see that the other two functions would initially appear as the better choice for minimizing error. Therefore, it is always important to decide what matters more to your program: stability over long time scales, energy conservation, and/or the number of steps it has to calculate to get a sufficiently small error.   

## Extensions

### Virial Theorem
TBA

## Languages, Libraries, Lessons Learned

The main library I learned and used in this project was the Scipy library. I found out about the different integration methods you can call from this library. Additonally, in the future I will have to consider the error accumulated from the method I use for simulation projects and whether the conservation of energy in the system would be an issue.

## Timekeeping

4/28: 22 hours

## Sources

Main Source: https://math.libretexts.org/Bookshelves/Differential_Equations/Numerically_Solving_Ordinary_Differential_Equations_(Brorson)/07%3A_Symplectic_integrators 

https://www.gorillasun.de/blog/euler-and-verlet-integration-for-particle-physics/ (Verlet and Symplectic Integrator)

https://en.wikipedia.org/wiki/Semi-implicit_Euler_method (Semi-implicit Euler wiki)

https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html (Modeled equation based of of def lotkavolterra from this section. I needed to pass t and return it as a  tuple for it to work.

ChatGPT (My Visual Studio stopped working at the beginning of the project and wouldn’t compile the python code or run it. I fed ChatGPT the error message and I eventually found a work around.) 

Prompt:

PS C:\Users\wolf1\OneDrive\Documents\GitHub\PHYS4130-S26> & C:/Users/wolf1/AppData/Local/Programs/Python/Python39/python.exe
Python 3.9.13 (tags/v3.9.13:6de2ca5, May 17 2022, 16:36:42) [MSC v.1929 64 bit (AMD64)] on win32
Type "help", "copyright", "credits" or "license" for more information.
Ctrl click to launch VS Code Native REPL
>>> plt.plot()
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
NameError: name 'plt' is not defined
>>> plt.plot()
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
NameError: name 'plt' is not defined
>>> v_array.append(start_v)
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
NameError: name 'v_array' is not defined
>>> import matplotlib.pyplot as plt
>>> 

https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.RK45.html (For RK45 solver’s notation)

https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.odeint.html (For Odient)

https://stackoverflow.com/questions/18648626/for-loop-with-two-variables (For loop zip)

https://en.wikipedia.org/wiki/Harmonic_oscillator (Harmonic Oscillator)
