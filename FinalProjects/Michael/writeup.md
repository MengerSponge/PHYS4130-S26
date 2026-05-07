## Numerical Methods for Time-Dependent PDEs
Partial differential equations (PDEs) arise naturally in modeling physical systems such as diffusion, wave propagation, and countless outher phenomena. While you can sometimes find analytical solutions, most realistic problems have nonlinear terms or complex boundary conditions that prevent frinding a closed form. Because of that, numerical methods are essential for approximating solutions and studying the behavior of these systems.

In this work, we examine two approaches for solving time-dependent PDEs. The Method of Lines converts the PDE into a system of ordinary differential equations (ODEs) by discretizing space, allowing the use of standard ODE integrators. We then introduce Exponential Time Differencing RK4 (ETDRK4), a more advanced method designed to handle stiff and higher-order problems efficiently. Together, these methods provide a foundation for numerically exploring a wide range of problems.

## The Method of Lines
### Theory
The first method we will investigate is known as the Method of Line. It combines the previously studied numeric differentiation stencils with ODE algorithms. Because of this, a problem must be written such that it is first order in time, and that time derivative should be isolated from the other terms of the PDE. In general for this method, you end up turning the PDE into a system of coupled ODEs in one of the coordinates of the system.

The process is best demonstrated by looking at an example. We will consider the heat equation in one dimension. 
```math
\partial_t u =\alpha^2 \, \partial_x^2 \, u 
```
The first step is to discretize our problem over our spatial coordinates. 
```math
\begin{gathered}
\text{ Let our spatial domain range from 0 to L. Let n be the number of discretized points.} \\
\text{ Thus, our domain points become} \quad D = \{0, \frac{L}{n}, \frac{2\,L}{n}, \dots , L\}
\end{gathered}
```
Then, we must discretize our solution over those points. 
```math
\begin{gathered}
\text{ Let } \, u_k(t) = u(x_k, t) \, \text{ where } \, x_k \in D
\end{gathered}
```
Now, we approximate the second derivative with a stencil using our spatial step size, h
```math
\begin{gathered}
h = L/n \\
\text{around} \, x_k, \, \partial_x^2 u \approx \frac{-u(x_k+2h, t) + 16\,u(x_k + h, t) - 30\,u(x_k, t) + 16\,u(x_k-h, t) - u(x_k-2h,t)}{12\,h^2}
\end{gathered}
```
Moving left or right by a step of size h is equivalent to moving left or right by an index in our discretized solution set. So, we obtain
```math
\begin{gathered}
\partial_x^2 u \approx \frac{-u_{k+2}(t) + 16\,u_{k+1}(t) - 30\,u_k(t) + 16\,u_{k-1}(t) - u_{k-2}(t),t)}{12\,h^2}.
\end{gathered}
```
Thus, we can put this result into the heat equation and obtain an equate the time derivative of each point in our solution to a function of the points around it. 
```math
\begin{gathered}
\frac{d}{dt}u_k(t) =  \frac{-u_{k+2}(t) + 16\,u_{k+1}(t) - 30\,u_k(t) + 16\,u_{k-1}(t) - u_{k-2}(t))}{12\,h^2}.
\end{gathered}
```
We are now done. Instead of a single function of two variables, we now have a system of functions in one variable that are all related through a coupled ODE system. At this point, we can churn this system of ODEs through a known integrator. Each discretized function serves as our approximation of the true solution around that point. The error for this method then comes from the order of your derivative stencil and the order of you ODE integrator.

This proces naturally generalizes to higher dimensions. Instead of a single index, we now have two or more (one for each spatial coordinate), and we have to approximate spatial derivatives in different directions. This method can also work for problems that are second order in time granted that you turn the problem into a system of PDEs that are first order in time. We demonstate this now. Consider a PDE such as 

```math
\begin{gathered}
\partial_t^2 \, = F(u,x,t). \\
\end{gathered}
```

We can transform this into a first order system through defining a new function, V.

```math
\begin{gathered}
\text{Let } \, v = \partial_t \, u \\
\text{Thus }\, \partial_t \begin{pmatrix} u \\ v \end{pmatrix} = \begin{pmatrix} v \\ F(u,x,t) \end{pmatrix}
\end{gathered}
```
Of course, now you have to solve for the other function over time as well by allowing a coupling of the two solutions. 

### Implementation
The first thing to decide is how the solutions will be stored. Numpy arrays are a natural choice. They allow the indexing of elements in the same way as the construction of the method. Furhtermore, this allows access to the vectorization of numpy objects, which will allow us to easily compute changes in our solution when time stepping. 

Next, some choices need to be made about the scope of the program. Since the solutions to the 2D versions of equations are more visually intersting than their one dimensional counterparts, and because we would like to be able to at least solve the wave equation, a reasonable level of generality to aim for is a program that can numerically solve a coupled system of two PDEs in two spatial dimensions that are first order in time. That is not a trivial task, and so we proceed carefully by first discretizing and initializaing our two solutions, U and V.
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
Then, to start time stepping we will need easy ways to compute various derivatives of these functions. Something that has not been addressed yet is the question of boundary conditions. For the time being, the program will be built presuming periodic boundary conditions (so inidicies loop back over to 0 when going beyond the array size), but the functionality will be implemeneted later on. Since we are using periodic BCs, we can use numpy's roll function to cleanly move around our arrays for computing the derivatives. For the equations we will be testing with this method, we only need the laplacian, but others can be implemented with some thought. We will use a nine point laplacian stencil for its 4th order error and its resistance to anisotropic effects. 
``` math
\begin{gathered}
\begin{pmatrix} 1/6 & 4/6 & 1/6 \\ 4/6 & -20/6 & 4/6  \\ 1/6 & 4/6 & 1/6\end{pmatrix}
\end{gathered}
```
The elements of this matrix represent the weights of the adjacent values in the sum to approximate the laplacian. The python implementation below utilizes the roll function to take in an array and computes the laplacian at all points.
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
Then, another useful function would be something that computes the derivatives for the two solutions. Here is an example for the wave equation.
```python
def F(t, U, V): #sample derivatives for the wave equation
    dU_dt = V
    dV_dt = Laplacian_9pt(U, h)

    return [dU_dt, dV_dt]
```
Now, all the structure is in place to compute things for this system. The only thing that remains is to package our system in a way the SciPy integrators can interpret. Currently, our system is represented as two N*N arrays whereas solve_ivp() needs the system to be a one dimensional array. We therfore will convert our data into that form. They can be rearranged with the flatten and concatenate functions.

```python
# Package the full system as one long array of all of our coupled solutions
Y0 = np.concatenate([U0.flatten(), V0.flatten()])
```

We will also need to compute the derivatives in this new representation. We will achieve this with some more array rearranging. We take in the packaged system, extract the correct forms for U and V, pass them through the derivative function, and then return the derivatives in the form we ne need for SciPy.

```python
def rhs(t, Y): #this computes the derivitives for the single list representation

    #First, extract U and V in their 2D array forms
    U = Y[:N*N].reshape((N, N)) #slice up to N*N - 1
    V = Y[N*N:].reshape((N, N)) #slice ater N*N
    
    # compute derivatives...
    dU_dt, dV_dt = F(t, U, V)

    #repackage our derivative to be in the same form factor as packaged system.
    return np.concatenate([dU_dt.flatten(), dV_dt.flatten()])
```

Now, we can start generating solutions by passing through solve_ivp(). We have some choices for our integrator. For these examples, we will go with RK45.

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
Where Nt is the number of frames we wanted to compute. Finally, we can start creating animated simulations for some PDEs. A good starting example is the heat equation.
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
In 2D, with periodic BCs, our simulation makes a solution that looks like this

https://github.com/user-attachments/assets/fb8d1d12-e28e-41f9-8485-cb393b5ffbc6

This looks very reasonable! The heat equation is the canoncial example of a diffusion system, which is exactly the behavior the simulation has captured. The next problem we can to try is the wave equation. 
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
So, V in our code represents the velocity at each point of U. This is expressed as follows.
```python
alpha = 2
def F(t, U, V):
    dU_dt = V
    dV_dt =  alpha**2 * Laplacian_9pt(U, h)
    return [dU_dt, dV_dt]
```
Then, with the same conditions as before, we get this animation. 

https://github.com/user-attachments/assets/5bda8ee5-9267-447b-8ea5-72c84b9ee01b

Which certainly looks like waves propogating. One last thing we can do to improve the program is add the ability to do fixed value boundary conditions. This can be done in the rhs function that is fed into solve_ivp. 
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
Where "boundary" is an array of points that is non-zero where we want to fix our values. The first round of updates is to ensure the values of the functions are fixed. The last round is to ensure that the time derivatives are zero to prevent the functions from attempting to update. If that were not there, then systems would slowly leak energy even if they theoreitcally shouldn't. By getting creative with how you set the boundaries, you can see how the wave equation behaves with a double slit and demonstrate that classical waves (somehow) behave as quantum particles. 

https://github.com/user-attachments/assets/3f4030f7-6b8d-415d-a6ce-1cd7f2b63f25

The last PDE we will examine is the 2D Kuramoto-Sivashinsky Equation (KSE). It's a popular example of a (relatively) simple PDE that exhibits spatiotemporal chaos. It's time evolution features the creation, intereaction, and annihiliation of small cell-like structures that bob around the simulation.

```math
\begin{gathered}
\partial_t \, u = - \nabla^2 u -  \nabla^4 u - |\nabla u|^2 \\
\text{where } \, \nabla^4 = \partial_x^4 + 2\, \partial_x^2 \, \partial_y^2 + \partial_y^4  \, \text{ is the bilaplacian}
\end{gathered}
```

Those 4th order derivatives and nonlinear terms are unsettling. Not to mention that it would be an ordeal to implement a function to compute the bilaplacian alone due to that mixed term. Even if we did do that, this equation is famously stiff. An algorithm as "basic" as the the method of lines can't possibly be good enough for this PDE with a resonable discretization since higher order stencis can be unstable. Therefore, we need a sophisticated method that can handle the numerous high-order spatial derivatives in a more precise manner than this while still having robust time stepping. Thankfully, such methods do exist, and we will see how one is constructed.

## Exponential Time Difference RK4 (ETDRK4)
### Theory
Let's briefly turn our attention from the KSE and instead look at the broader class of problems that it falls under. Condsider PDEs of the form

```math
\begin{gathered}
\partial_t u = Lu + N(u,t) \\
\text{where L is a linear operator and N(u,t) is the nonlinear term} \,
\end{gathered}
```

subjected to periodic boundary conditions. To simulate this problem, we will construct an algorithm known as Exponential Time Difference RK4 (ETDRK4). It is a pseudo spectral method that is designed for handling stiff, high order problems. The derivation starts with creating an integrating factor.

```math
w = e^{-Lt} \,u.
```

Recall that the operator exponential is itself an operator. We seek an update formula for u. We do so by first differentiating w.

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

Where in that last step we have made a change of varaibles. Now, we substitute in our definition of w.

```math
\begin{align*}
e^{-L\,(t+h)} \,u - e^{-Lt} \,u = \int_{0}^{h} e^{-L\,(t+\tau)} \, N(u(t+\tau),t+\tau) d\tau \\
\end{align*}
```

By shuffling around terms and canceling out repeating factors, we obtain an exact formula for the change in u over a time step of size h.

```math
u(x,y, t+h) = e^{Lh}\,u(x,y,t) + e^{Lh}\, \int_{0}^{h} e^{-L\tau}\, N(u(x,y,t+\tau,t+\tau)\,d\tau.
```

By creating some discretization variables, we can rewrite this as an update formula.

```math
\begin{align*}
&\text{Let} \, h ,\ \text{be a time step} \\
&\text{Let} \, t_n \, = n\,h \\
&\text{Let} \, u_{n+1} = e^{L\,h}\,u_n + e^{L\,h}\, \int_{0}^{h} e^{-L\tau}\, N(u(x,y,t+\tau,t+\tau)\,d\tau.
\end{a
lign*}
```

Before proceeding, we should note some that we never discetized the operator L. In fact, if there is no nonlinear term, then this method will be exact. We only need to approxmiate in the nonlinear step, which we adress next. 

Up until this point we have not done anything different from the construction of other exponential time differencing methods. The distinction comes in how you approximate the integral of the nonlinear term. Hence the name, ETDRK4 uses an RK4-like quadrature to approximate it in terms of our discretization variables. The derivation is not presented here, but can be found in the original paper detailing this method titled "Exponential Time Differencing for Stiff Systems" by Cox and Matthews. Our quadtrature coefficients are computed with the following formulas.

```math
\begin{align*}
&f_1 ​= (L\,h)^{−2} \, (− 4 − L\,h + e^{L\,h}\,(4−3\,L\,h+(L\,h)^2)) \\
&f_2 = (L\,h)^{−2} \, (2 + L\,h + e^{L\,h}\,(−2+ L\,h)) \\
&f_3 ​= (L\,h)^{−2}(− 4 − 3\,L\,h − (L\,h)^2 + e^{L\,h}\,(4 − L\,h)) \\ \\

&a = e^{L\,h/2}\,u_n + L^{-1} \,(e^{L\,h/2} - I)\,N(u_n) \\
&b = e^{L\,h/2}\,u_n + L^{-1} \,(e^{L\,h/2} - I)\,N(a) \\
&c = e^{L\,h/2}\,a  +  L^{-1} \,(e^{L\,h/2} - I)\,(2\,N(b) - N(u_n)) \\
\end{align*}
```

We put these together for our final update formula. 

```math
\begin{align*}
&u_{n+1} = e^{L\,h}\,u_n + h\,(N(u_n)\,f_1 + 2\,(N(a) + N(b)) \, f_2 + N(c) \, f_3) \\
\end{align*}
```

At this point, the theory for this method has technically been fully developed. You can, hypothetically, go and compute all the update terms by expanding these operator exponentials. This is, of course, silly to do in practice. One goal for deriving a numerical method is to have it be computationally efficient. Therefore, we develop a little more theory that helps in computing these coefficients. We proceed by writing our function as a Fourier series over 2D plane waves. (we may do this because plane waves satisfy periodic BCs)

```math
u(x,y,t) = \sum_{k_x , k_y} A_{k_x, k_y}(t) e^{i(k_x  x + k_y  y)},
```

Then, because L is a linear differential operator, it immediatley follows that it is diagonal in this Fourier space representation of U. This is because the plane waves are the eigenstates of differentiation. As an example, consider the partial derivative with repect to x on a generic Fourier mode. 

```math
\begin{align*}
\partial_x \, e^{i\,(k_x \, x + k_y y)} &= i\,k_x\,e^{i\,(k_x \, x + k_y y)}
\end{align*}
```

This will become useful when we impliment the algorithm into python. Thanks to the Fast Fourier Transform (FFT), we can efficiently go from the position space to the frequeny space representations of the function. Therefore, our computations of operator exponentials reduce from demanding matrix multiplication to simple evaluations of the exponential at several points. Again, we demonstrate this with the partial derivative with respect to x. 

```math
\begin{align*}
e^{\partial_x} \, e^{i\,(k_x \, x + k_y y)}
&= \sum_{n = 0}^{\infty} \frac{\partial_x ^n }{n\,!} \, e^{i\,(k_x \, x + k_y y)} \\
&= \sum_{n = 0}^{\infty} \frac{i^n \, k_x^n }{n\,!} \, e^{i\,(k_x \, x + k_y y)} \\
&= e^{i \, k_x} e^{i\,(k_x \, x + k_y y)} \\
\end{align*}
```

Therefore, any time we need to compute something involving a derivative we can simply do it in Fourier space. The very final piece of theory that needs to be built has to do with the computaion of the coefficients f1, f2, and f3 for the update formula. Those will not vary as we step through time. When we evalutate them in fourier space, we just substitute in the eigenvalues of L as we did for differentiation. However, these eigenvalues can be in regions where the formulas risk catastrophic cancellation or divergence. Therefore, they can be computationally unstable. Thankfully,a trick devleloped by Kassam and Trefethen addresses this. Recall Cauchy's Formula for integration from Complex Analysis. 

```math
f(a) = \frac{1}{2\pi i} \, \oint_{\gamma} \frac{f(z)}{z - a}\, dz
```

We will use this to evaluate f1, f2, and f3. Since our precarious points lie near 0, we take our contour to be the unit circle centered around our evaluation point, a. Then, we can just use trapezoidal rule to approximate the integral. 

```math
f(a) \, \approx \, \frac{2\pi}{N} \sum_{k=0}^{N-1} f\!\left(a + e^{2\pi i k / N}\right)\, i e^{2\pi i k / N}
```

In practice, this may leave small imaginary terms in our computation. We will take only the real value of this since we are only interested in real-valued problems for this. Now that all the theory is complete, we can see how this is implemented in python. Surprsingly, despite being more theoreitcally complex, the implimenttion will be considerably simpler than the method of lines. 

### Implimentation
Since we were interested in the KSE, we will implement our program with that problem in mind. We can indentify our L and N terms for ETDRK4 as 

```math
\begin{aligned}
L\, = \, - \nabla^2 u - \nabla^4 u \\
N(u) = - |\nabla u|^2 \\
\end{aligned}
```

Furthermore, it is easy to show that 

```math
L\,e^{i\,(k_x \, x + k_y y)} \, = \, (k_x^2 + k_y^2 - (k_x^2 + k_y^2)^2)\,e^{i\,(k_x \, x + k_y y)}
```

With those two things in mind, we begin putting together the program. We start by initializing our postion, frequency, and L arrays as well as initial conditions on our solution.

```python
# N = number of points along a spatial dimension
N = 100

# size = length along the spatial dimensions
size = 75

# spatial arrays (uncombined) will be combined into an array of coordinates
x = size*np.arange(0,N)/N
y = size*np.arange(0,N)/N

#our space
X, Y = np.meshgrid(x, y, indexing='ij')

dx = x[1] - x[0]
dy = y[1] - y[0]

kx = np.fft.fftfreq(N, d=dx) * 2*np.pi
ky = np.fft.fftfreq(N, d=dy) * 2*np.pi

# our array of eigenvlaues of L
Kx, Ky = np.meshgrid(kx, ky, indexing='ij')

K2 = Kx**2 + Ky**2

#fourier space matrix for L
L = 0.25*(K2 - K2**2) 

# random noise IC
u0 = 0.01*np.random.randn(len(x), len(y))
```

Then, we set up code for things we can precompute.

```python
t_steps = int(Tmax/h)
# Operator Exponentials 
E = np.exp(h*L)
E2 = np.exp(0.5*h*L)

# Set up to do contour integral trick
M = 20
r = np.exp(1j*np.pi*(np.arange(1, M+1) - 0.5)/M) #pieces of the unit circle
z = h * (L) #arguments of the phi functions
LR = z[:, :, None] + r[None, None, :] #take all the points in z and add a unit circle around them

# the actual coefficents computed using the contour integral trick
f1 = h*np.mean(((-4-LR+np.exp(LR)*(4-3*LR+LR**2))/LR**3),axis=2).real
f2 = h*np.mean(((2+LR+np.exp(LR)*(-2+LR))/LR**3),axis=2).real
f3 = h*np.mean(((-4-3*LR-LR**2+np.exp(LR)*(4-LR))/LR**3),axis=2).real

Q  = h*np.mean(((np.exp(LR/2) - 1) / (LR)), axis=2).real
```

Next, we need to set up a function that takes in our Fourier space representation of our function, computes the nonlinear term in positon spsace since that is harder, and then returns the Fourier space representation of the nonlinear term.

```python
def NonLinear(u_ft, t): #takes in fourier transformed u and computes the non linear part.
    #spatial derivatives
    u_x = np.fft.ifft2(-1j*Kx*u_ft).real
    u_y = np.fft.ifft2(-1j*Ky*u_ft).real

    # nonlinear terms in position space
    N = -0.25 * (u_x**2 + u_y**2)
    return np.fft.fft2(N) #return to fourier space
```

Note here that we are using the 2D FFT function in numpy because our problem is 2D. Now, we can define a function to do a single time step in fourier space. 

```python
def Step(u_ft, t):
    t = n*h

    #RK4 stepping
    Nu = NonLinear(u_ft, t)
    a = E2*u_ft + Q*Nu
    Na = NonLinear(a ,t + h/2)
    b = E2*u_ft + Q*Na
    Nb = NonLinear(b, t + h/2)
    c = E2*a + Q*(2*Nb-Nu)
    Nc = NonLinear(c ,t + h)
    u_ft_new = E*u_ft + Nu*f1 + 2* (Na+Nb)*f2 + Nc*f3

    return u_ft_new
```

Finally, we can easily create a loop to compute our solution. 

```python

#our initial state in fourier space
v = np.fft.fft2(u0) #create v as the fourier transform of u0

#start these lists to save the funciton over time
u_list = [u0]
t_list = [0]

#nplt = number of steps stored in the animation. nplt = 1 means save all. larger means save less
nplt = 13

#Tmax = maximum time
Tmax = 300

# h = time step size
h = 0.25

#Verbosity controls how often the program reports progress
verbosity = int(Tmax/(10*h))

for n in range(1,t_steps+1):
    #step
    t = n*h
    v = Step(v, t)

    if n%nplt == 0: #to save for animating
        u = np.fft.ifft2(v).real
        u_list.append(np.copy(u))
        t_list.append(np.copy(t))
    if n%verbosity == 0: # to output progress
        print(100*n/len(range(1,t_steps+1)), " % complete")
````

When we test this out with with the KSE, our animated solution looks like this.

https://github.com/user-attachments/assets/2c8d3b6c-598c-4257-ae0e-1bc2f842e83f

We can kind of see our expected dynamics, but there is this large uniform background that divereges as time goes on. To trouble shoot this, we need to look at any terms in our solution that correspond to a constant term in space. In the fourier series, we see that this term goes with the (0,0) mode. Then, let's compute how this constant term evolves according to the KSE. 

```math
\begin{aligned}
L \, A_{(0,0)} (t) \, e^{i\,(0*x + 0*y)} = L \, A_{(0,0)} (t) = 0 \, \text{ since L is a spatial differential opeartor} \\
N(A_{(0,0)} (t) \, e^{i\,(0*x + 0*y)}) = N(A_{(0,0)} (t)) = 0 \, \text{ since N(u) is composed of spatial differential opeartors} \\
\text{Thus, the (0,0) mode is constant in time. }
\end{aligned}
```

Since our divergence is uniform over space, it seems that the simulation has the wrong time evolution for the constant term. Then, we also observe that our solution diverges slowly over the simulation time. Therefore, we conclude that our approximation for the dynamics of the (0,0) mode is just slightly off form zero in the same direciton at each step. Thus, we accumulate these little errors over enough time to see diveregence. This is an occasional problem with tihs method, and it is easily fixed by enforcing the (0,0) mode to be constant in time inside our step function. 

```python
def Step(u_ft, t):
    a00 = u_ft[0,0] #to enforce constant average
    t = n*h

    #RK4 stepping
    Nu = NonLinear(u_ft, t)
    a = E2*u_ft + Q*Nu
    Na = NonLinear(a ,t + h/2)
    b = E2*u_ft + Q*Na
    Nb = NonLinear(b, t + h/2)
    c = E2*a + Q*(2*Nb-Nu)
    Nc = NonLinear(c ,t + h)
    u_ft_new = E*u_ft + Nu*f1 + 2*(Na+Nb)*f2 + Nc*f3

    #enforce constant average
    u_ft_new[0,0] = a00
    return u_ft_new
```

Now, with this fix implemented, our simulations are exactly creating the solutions that we expect. 

https://github.com/user-attachments/assets/628cd554-8723-4b8c-90b0-6d4a254f8e1d

At last, we have achieved our goal. As one last reward for slogging through the math, here are two KSE simulations with highly symmetric ICS.

https://github.com/user-attachments/assets/c0f20e9d-e0f4-4025-aab6-243ab9ce39be

https://github.com/user-attachments/assets/552360a9-2d23-4142-abf9-c47f59d31d41

We see that, for a chaotic system, this method preserves the symmetry for quite a long time! Of course, because error is unavoidable, eventually the symmetry breaks down and the system devolves into asymmetric chaos. 

## Conclusion 
Two numerical methods for solving PDEs have been demonstrated and applied to some example problems. The first one, the Method of Lines, use numeric differentiation stencils to approximate the PDE as a coupled system of ODEs that can then be solved numerically using any known integrator. This method was implemented using python and shown to work for two simple PDEs in the heat and wave equation. Then, the derivation was outlined for Exponential Time Difference RK4. It was implemented and tested with the KSE, demonstrating the effectiveness of this method for high-order and stiff problems. With both of thse methods, many problems can be explored.

## Time Keeping
I spent about 30 hours on the program and 15 hours on the report.

## Languages, Libraries, Lessons Learned
The only language used was python. Specific libraries were numpy, scipy.integrate ,matplotlid.pyplot, and matplotlib.animation. From numpy, arrays and the FFT function was used. scipy.integrate was used in the method of lines for solving the ODE system, and matplotlib was used for plotting and animating. The main skill I developed from this, beyond learning to numerically solve PDEs, was how to be clever when indexing. 

## Attribution
The following works, while not directly cited, were essential in putting together my own derivation outline of ETD RK4 (I already knew of the method of lines from personal projects.)

1) "Exponential Time Differencing for Stiff Systems" by S. M. Cox and P.C. Matthews
2) "Hydrodynamics of the Kuramoto Sivashinsky Equation in Two Dimensions" by B. M. Boghosian et al.
3) "An Exponential Time Differencing Scheme for the Kuramoto-Sivashinksy Equation" by G. Zavalani

Also of note is the instagram account "MATH_INFINITUM" for the reel they posted on the Kuramoto-Sivashinksy Equation.
