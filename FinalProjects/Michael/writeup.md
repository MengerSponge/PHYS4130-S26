## Numerical Methods for Time-Dependent PDEs
dummy intro
There are PDEs that must be solved. But, they don't want to be solved. However, mathematicians are clever. THey found a way to solve problems that don't want to be solved!

## The Method of Lines: Theory
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

## The Method of Lines: Implementation
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
Then, to start time stepping we will need easy ways to compute various derivatives of these functions. The basic ones are the first order partial derivatives and the laplacian. Something that has not been addressed yet is the question of boundary conditions. For the time being, the program will be built presuming periodic boundary conditions (so inidicies loop back over to 0 when computing coupled terms), but the functionality will be implemeneted later on. Since we are using periodic BCs, we can use numpy's roll function to move around in our arrays for computing the coupled derivatives. For the simulations we will be testing with this method, we only need the laplacian, but others can be implemented with some thought. We will use a nine point laplacian stencil for its higher error order (4th order) and its resistance to anisotropic effects. 

``` math
\begin{gathered}
\begin{pmatrix} 0.25 & 0.5 & 0.5 \\ 0.5 & -3 & 0.5  \\ 0.25 & 0.5 & 0.5\end{pmatrix}
\end{gathered}
```
Where the elements of this matrix represent the weights of the adjacent values in the sum to approximate the laplacian
```python
def Laplacian_9pt_2(U, h): #different stencil

    return (
        2*(np.roll(U, 1, 0) + np.roll(U, -1, 0) + np.roll(U, 1, 1) + np.roll(U, -1, 1))
        +
        (
            np.roll(np.roll(U, 1, 0), 1, 1)+
            np.roll(np.roll(U, 1, 0), -1, 1)+
            np.roll(np.roll(U, -1, 0), 1, 1)+
            np.roll(np.roll(U, -1, 0), -1, 1)
        )
        - 3*U
    )/(2* h**2)
```







## Old Stuff from the outline
a) There is a hihgly intuitive way of numerically solving PDEs. 
  I) We are able to discretize spatial derivativres over our spatial domain and create a 
     large system of coupled ODEs since the approximation of the spatial derivatives
     mixes discretized terms.
  II) Then, we can churn this system of ODEs through a numeerical method such as RK45 
      and create a numerical solution. 
  III) Here, show a plot of a sample simulation that I can animate
b) That method is fine enough for simple PDEs, but if you want to do something 
   more complicated then you need a more sophisticated method. If we want to solve something
   like the Kuramoto-Sivashinksy equation (which is 4th order in space and nonlinear) then the old method
   won't cut it. Instead, we can use Exponential Time Difference RK4 (ETDRK4)
   I) Explaining the method: We want to solve a PDE of the form 
```math
\partial_t u = Lu + N(u)
```
  Where L is a linear (but possibly stiff) operator and N(u) is nonlinear. 

  II) Then, if L is a differential operator, and we solve our problem over periodic BCs, then we can write our funciton as a fourier expansion
```math
u(x,y,t) = \sum_{k_x , k_y} A_{k_x, k_y}(t) e^{i(k_x  x + k_y  y)},
```
  and L is now a diagonal operator in the sense that u is represented in the eigenbasis of L. 
  Let 
```math
\hat{u(t)} = \begin{pmatrix} A_{k_x, k_y}(t) \end{pmatrix}_{k_x, k_y wave numbers}
```
  be the fourier space representation of u. 
  III) Let h be a timestep. We now introduce an integrating factor
```math
w(x,y,t) = e^{-Lh}.
```
  Recall that the matrix exponential is itself an operator. This allows us to derive an exact time step formula.
  (Skipping steps here to get to the time step formula)
```math
u(x,y, t+h) = e^{Lh}u(x,y,t) + e^{Lh} \int_{0}^{h} e^{-L\tau} N(u(x,y,t+\tau)d\tau.
```
This formula is exact, and the RK4 part comes in how we approximate this integral (not included here atm)
We need the fourier space representation for the efficient computation of matrix exponentials. 

# Solving the Kuramoto Sivashinksy Equation
Derive the particular update formula for ETDRK4 here for the KS equation. Show pretty plots, and the python implementation here.
  
  
