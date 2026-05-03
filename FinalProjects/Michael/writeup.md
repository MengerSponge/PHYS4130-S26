## Numerical Methods for Time-Dependent PDEs
dummy intro
There are PDEs that must be solved. But, they don't want to be solved. However, mathematicians are clever. THey found a way to solve problems that don't want to be solved!

## The Method of Lines
### Theory
The first method we will investigate is known as the method of lines. It builds off of already devloped ODE  integrators as well as the numeric differentiation stencils that we studied in previous projcets. Because of this, the problem must be written in a form such that it is first order in time, and that time derivative should be isolated from the other terms of the PDE. The goal is to turn the PDE into a system of coupled ODEs in all but one coordinate. 

The process is best demonstrated by first looking at an example. We will consider the heat equation in one dimension. 
```math
\partial_t u = \partial_x^2 \, u 
```
The first step is to discretize our problem over our spatial coordinates.
```math
\begin{gathered}
\text{ Let our spatial domain range from 0 to L. Let n be the number of discretized points.} \\
\text{ Thus, our domain points become} \quad D = \{0, \frac{L}{n}, \frac{2\,L}{n}, \dots , L\}
\end{gathered}
```
Then, we must discretize our solution over our points. 
```math
\begin{gathered}
\text{ Let } \, u_k(t) = u(x_k, t) \, \text{ where } \, x_k \in D
\end{gathered}
```
Now, we need to approximate the second derivative with a stencil using our step size 
```math
\begin{gathered}
h = L/n \\
\text{around} \, x_k, \, \partial_x^2 u \approx \frac{-u(x_k+2h, t) + 16\,u(x_k + h, t) - 30\,u(x_k, t) + 16\,u(x_k-h, t) - u(x_k-2h,t)}{12\,h^2}
\end{gathered}
```
Moving left or right by a step of size h is equivalent to moving left or right by an index in our discretized solution. So, we obtain
```math
\begin{gathered}
\partial_x^2 u \approx \frac{-u_{k+2}(t) + 16\,u_{k+1}(t) - 30\,u_k(t) + 16\,u_{k-1}(t) - u_{k-2}(t),t)}{12\,h^2}.
\end{gathered}
```
Thus, we can put this result into the heat equation and obtain an equation for the time derivative of each point in our solution as a funciton of the points around it. 
```math
\begin{gathered}
\frac{d}{dt}u_k(t) =  \frac{-u_k+2(t) + 16\,u_k+1(t) - 30\,u_k(t) + 16\,u_k-1(t) - u_k-2h(t),t)}{12\,h^2}.
\end{gathered}
```
We are now done. Instead of a single function of two variables, we now have a system of functions in one vartiable that are all related through a coupled ODE system. At this point, we then churn this system of ODEs through a known integrator. Each discretized function serves as our approximation of the function around that point in space. This method is modular in the sense that you select the order of your spatial error by choosing the order of your derivative stencil. 

This proces naturally generalizes to higher dimensions. Instead of a single index, we now have two or more (one for each spatial coordinate), and we have to approximate spatial derivatives in the correct direction. This method can also work for problems that are second order in time granted that you turn the problem into a system of PDEs that are first order in time, which will be shown below. Consider a PDE such as 
```math
\begin{gathered}
\partial_t^2 \, = F(u,x,t). \\
\end{gathered}
```
We can transform this into a first ordey system through the following
```math
\begin{gathered}
\text{Let } \, v = \partial_t \, u \\
\text{Thus }\, \partial_t \begin{pmatrix} u \\ v \end{pmatrix} = \begin{pmatrix} v \\ F(u,x,t) \end{pmatrix}
\end{gathered}
```
Of course, now you have to solve for the other function over time as well by allowing a coupling of the two through derivatives. This also allows you to consider systems of PDEs that are first order in time. 

### Implementation
The first thing to decide is how the solutions will be stored. Numpy arrays are a natural choice. They allow the indexing elements in the same way as the construction of our method. Furhtermore, this allows acces to the vectorization of numpy arrays, which will allows us to easily compute changes in our solution when time stepping. 

Then, some choices need to be made about the scope of the program. Since the 2D versions of famous equations such as the heat or wave equation have more visually intersting solutions than their one dimensional counterparts, and because we would like to be able to at least solve the wave equation, a reasonable level of generality to aim for is a program that can numerically solve a coupled system of two PDEs in two spatial dimensions that are first order in time. That is not a trivial task,and so we proceed carefully by first discretizing and initializaing our two solutions: U and V.
```python

N = 151 # Num points on a postion side
h = 0.5 # spatial step size

T_final = 500
Nt = 800

dt = T_final/Nt #time step

# Sample Initial conditions
U0 = np.zeros((N, N))
U0[120:130, N//2 - 5:N//2 ] = 20

V0 = np.zeros((N, N))
```
Then, to start time stepping we will need easy ways to compute various derivatives of these functions. The basic ones are the first order partial derivatives and the laplacian. Something that has not been addressed yet is the question of boundary conditions. For the time being, the program will be built presuming periodic boundary conditions (so inidicies loop back over to 0 when computing coupled terms), but the functionality will be implemeneted later on. Since we are using periodic BCs, we can use numpy's roll function to move around in our arrays for computing the coupled derivatives. For the simulations we will be testing with this method, we only need the laplacian, but others can be implemented with some thought. We will use a nine point laplacian stencil for its 4th order error and its resistance to anisotropic effects. 
``` math
\begin{gathered}
\begin{pmatrix} 1/6 & 4/6 & 1/6 \\ 4/6 & -20/6 & 4/6  \\ 1/6 & 4/6 & 1/6\end{pmatrix}
\end{gathered}
```
Where the elements of this matrix represent the weights of the adjacent values in the sum to approximate the laplacian. The python implementation below utilizes the roll function to take in an array and computes the laplacian at all points.
```python
def Laplacian_9pt(U, h): #my old one (this) is faster

    return (
        4*(np.roll(U, 1, 0) + np.roll(U, -1, 0) + np.roll(U, 1, 1) + np.roll(U, -1, 1))
        +
        (
            np.roll(np.roll(U, 1, 0), 1, 1)+
            np.roll(np.roll(U, 1, 0), -1, 1)+
            np.roll(np.roll(U, -1, 0), 1, 1)+
            np.roll(np.roll(U, -1, 0), -1, 1)
        )
        - 20*U
    )/(6 * h**2)
```
Then, another useful function would be something that computes the derivatices for the two solutions.
```python
def F(t, U, V): #sample derivatives for the wave equation
    dU_dt = V
    dV_dt = Laplacian_9pt(U, h)

    return [dU_dt, dV_dt]
```
Now, all the structure is in place to compute things for this system. The only thing that remains is to find a way to package our system in a way the SciPy integrators can interpret. Currently, our system is represented as two N*N arrays whereas solve_ivp() needs the system to be a one dimensional array. With that in mind, we need to convert our data into that form. We will also need to compute the derivatives in this new representations accordingly. We will achieve this by flattening and concatenating our two arrays. 
```python
# Package the full system as one long array of all of our coupled solutions
Y0 = np.concatenate([U0.flatten(), V0.flatten()]) 

def rhs(t, Y): #this computes the derivitives for the single list representation

    #First, extract U and V in there 2D array forms
    U = Y[:N*N].reshape((N, N)) #slice up to N*N - 1
    V = Y[N*N:].reshape((N, N)) #slice ater N*N
    
    # compute derivatives...
    dU_dt, dV_dt = F(t, U, V)

    #package our derivative to be in the same form factor as packaged system.
    return np.concatenate([dU_dt.flatten(), dV_dt.flatten()])
```
Now, we can start generating solutions by passing through solve_ivp(). We have some choices for our integrator. For this, we will go with RK45.
```python
from scipy.integrate import solve_ivp

sol = solve_ivp(
    rhs, 
    t_span = (0, T_final), 
    y0 = Y0,
    method = 'RK45',
    t_eval = np.linspace(0, T_final, Nt)
    )
```
We can then extract our solutions and reshape them accordingly to minimc our original structure. 
```python
U_list = [sol.y[:N*N,k].reshape((N,N)) for k in range(0, Nt)]
V_list = [sol.y[N*N:,k].reshape((N,N)) for k in range(0, Nt)]
```
Where Nt is the number of frames we wanted to compite. Finally, we can start creating animated simulations for some PDEs. A good starting example is the heat equation that was used when demonstrating the method. 
```math
\begin{aligned}
\partial_t \, u = \alpha^2 \, \nabla^2 \, u
\end{aligned}
```
We only need one function, so in our code we can simply update U while setting the derivatives for V to be an array of zeros. 
```python
alpha = 2
def F(t, U, V):
    dU_dt = alpha**2 * Laplacian_9pt(U, h)
    dV_dt = np.zeros_like(V)
    return [dU_dt, dV_dt]
```
In 2D, with periodic BCs and ICs fixing some non-zero temperature at certain point, we get a simulation like the following. 

![Heat Equation](Heat_Eqn_Periodic.gif)

This looks very reasonable! The heat equation is the canoncial example of a diffusion system, which is exactly the behavior the simulation has demonstrated. The next problem we can to try is the wave equation. 
```math
\begin{gathered}
\partial_t^2 \, u = \alpha^2 \, \nabla^2 \, u
\end{gathered}
```
We need to do some more work to adapt this to the program. The system must be expressed as a coupled set of PDEs that are first order in time. 
```math
\begin{gathered}
\text{Let } \, v = \partial_t \, u \\
\text{Thus }\, \partial_t \begin{pmatrix} u \\ v \end{pmatrix} = \begin{pmatrix} v \\  \nabla^2 \, u  \end{pmatrix}
\end{gathered}
```
So, V in our code represents the velocity at each point of U. The system is expressed in the following way.
```python
alpha = 2
def F(t, U, V):
    dU_dt = V
    dV_dt =  alpha**2 * Laplacian_9pt(U, h)
    return [dU_dt, dV_dt]
```
Then, with the same conditions as before, the simulation produces the following animation. 

![Wave Equation](Wave_Eqn_Periodic.gif)

Which definitley looks like waves propogating. One last thing we can do to improve the program is add the ability to do fixed value boundary conditions. This can be done in the rhs function that is fed into solve_ivp. The modified function is as follows.
```python
Do_BCs = True

def rhs(t, Y):
    
    U = Y[:N*N].reshape((N, N)) #slice up to N*N - 1
    V = Y[N*N:].reshape((N, N)) #slice ater N*N
    
    #Enforce values
    if Do_BCS == True:
        U[boundary]   = 0
        V[boundary] = 0
    
    # compute derivatives...
    dU_dt, dV_dt = F(t, U, V)

    #Enforce dynamics to stop them from updating
    if Do_BCS == True:
        dU_dt[boundary]   = 0
        dV_dt[boundary] = 0

    return np.concatenate([dU_dt.flatten(), dV_dt.flatten()])
```
Where "boundary" is an array of points that is non-zero where we want to fix our values. The first round of updates is to ensure the values of the functions are fixed. The last round is to ensure that the time derivatives are zero to prevent the functions from attempting to update. If that were not there then systems would slowly leak energy even if they theoreitcally shouldn't. By getting creative with how you set the boundaries, you can see how the wave equation behaves with a double slit and demonstrate that classical waves behave as quantum particles. 

![Wave Equation](Wave_Eqn_Boundary.gif)

Another interesting equation to simulate is the 2D Kuramoto-Sivashinsky Equation (KSE). It's a popular example of a simple PDE that exhibits spatiotemporal chaos. It's time evolution features the creation, intereaction, and annihiliation of small cell-like structures that bob around before fading away.
```math
\begin{gathered}
\partial_t \, u = \nabla^2 u + \nabla^4 u + |\nabla u|^2 \\
\text{where } \, \nabla^4 = \partial_x^4 + 2\, \partial_x^2 \, \partial_y^2 + \partial_y^4  \, \text{ is the bilaplacian}
\end{gathered}
```
Those 4th order derivatives and nonlinear terms are unsettling. Not to mention that it would be an ordeal to implement a function to compute the bilaplacian alone due to that mixed term. Even if we did do that, this equation is famously stiff. An algorithm as basic as the the method of lines can't possibly be accurate enough for this PDE with a resonable spatial or temporal discretization. Therefore, we need a sophisticated method that can handle the numerous spatial derivatives in a more precise manner than this while still having robust time stepping. Thankfully, such methods do exist, and we will see how one is constructed.

## Exponential Time Difference RK4 (ETDRK4)
### Theory
Let's briefly turn our attention from the KSE and instead look at the broader class of problems that it falls under. Condsider PDEs of the form
```math
\begin{gathered}
\partial_t u = Lu + N(u,t) \\
\text{where L is a linear operator and N(u,t) is the nonlinear term} \,
\end{gathered}
```
subjected to periodic boundary conditions. We now introduce an integrating factor
```math
w = e^{-Lt} \,u.
```
Recall that the matrix exponential is itself an operator. We seek to find an update formula for w that then can be turned into an update formula for u. Differentiating w yields
```math
\begin{align*}
\partial_t \, w &= \partial_t(\, e^{-Lt} \,u.) \\
                &= \partial_t(\, e^{-Lt})  \, u + e^{-Lt} \, \partial_t(\,u) \\
                &= -L\, e^{-Lt} \, u + e^{-Lt} \, (Lu + N(u,t)) \\
                &=  e^{-Lt} \, N(u,t).
\end{align*}
```
Then, we let h be a time step. Our exact change in w as we go from t to t+h comes by integrating.
```math
\begin{align*}
w(t+h) - w(t) &= \int_{t}^{t+h} \partial_t w dt \\
              &= \int_{t}^{t+h} e^{-Lt} \, N(u,t) dt \\
              &= \int_{0}^{h} e^{-L\,(t+\tau)} \, N(u(t+\tau),t+\tau) d\tau \\
\end{align*}
```
Where in that last step we have made a change of varaibles. Now, we substitute in our definition of u.
```math
\begin{align*}
e^{-L\,(t+h)} \,u - e^{-Lt} \,u = \int_{0}^{h} e^{-L\,(t+\tau)} \, N(u(t+\tau),t+\tau) d\tau \\
\end{align*}
```
By shuffling around terms and canceling out repeating factors, we obtain an exact formula for the change in u over a time step of size h.
```math
u(x,y, t+h) = e^{Lh}\,u(x,y,t) + e^{Lh}\, \int_{0}^{h} e^{-L\tau}\, N(u(x,y,t+\tau,t+\tau)\,d\tau.
```
Up until this point, we have not done anything different from the construction of other exponential time differencing methods. The distinction comes in how you approximate the integral. Hence the name, ETDRK4, we will use a RK4-like quadrature to approximate it. 
```math
u(x,y,t) = \sum_{k_x , k_y} A_{k_x, k_y}(t) e^{i(k_x  x + k_y  y)},
```
Let 
```math
\hat{u(t)} = \begin{pmatrix} A_{k_x, k_y}(t) \end{pmatrix}_{k_x, k_y wave numbers}
```
be the fourier space representation of u. 

# Solving the Kuramoto Sivashinksy Equation
Derive the particular update formula for ETDRK4 here for the KS equation. Show pretty plots, and the python implementation here.
  
  
