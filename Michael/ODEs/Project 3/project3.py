import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp, odeint

#Numerical Solvers
#As a convenience, the Euler and RK solvers will be configured to solve 
#dx/dt = f(t, x,y), dy/dt = g(t, x,y)
#Where x and y may be numpy arrays. This works because numpy functions are vectorized.
#Note: RK2,RK4, and Verlet algorithms soruced from: "Computational Physics" by Prof. Mark Newman
def Euler(X0,Y0, tmin, tmax, nts, f,g):
    t = np.linspace(tmin,tmax,nts,endpoint = False)
    X = np.zeros(nts)
    Y = np.zeros(nts)

    dt = t[1]-t[0]

    X[0] = X0
    Y[0] = Y0

    for it in range(0,nts):
        #Euler coefficients
        k1x = dt*f(t[it], X[it], Y[it])
        k1y = dt*g(t[it], X[it], Y[it])

        #Euler update
        X[it+1] = X[it] + k1x
        Y[it+1] = Y[it] + k1y
        
    return t, X, Y

def RK2(X0,Y0,tmin,tmax,nts,f,g):
    t = np.linspace(tmin,tmax,nts,endpoint = False) #time points
    X = np.zeros(nts)
    Y = np.zeros(nts)

    dt = t[1]-t[0]

    X[0] = X0
    Y[0] = Y0

    for it in range(0,nts):
        #RK2 Coefficients
        k1x = dt*f(t[it], X[it], Y[it])
        k1y = dt*g(t[it], X[it], Y[it])

        k2x = dt*f(t[it] + 0.5*dt, X[it] + 0.5*k1x, Y[it] + 0.5*k1y)
        k2y = dt*g(t[it] + 0.5*dt, X[it] + 0.5*k1x, Y[it] + 0.5*k1y)
        #RK2 update
        X[it+1] = X[it] + k2x
        Y[it+1] = Y[it] + k2y

    return t, X, Y 

def RK45(X0,Y0,tmin,tmax,nts,du_dt):
    t_span = (tmin,tmax)
    t = np.linspace(tmin,tmax,nts,endpoint = False)
    U0 = [X0,Y0]

    return solve_ivp(du_dt, t_span, U0,t_eval = t, method = 'RK45')

def LSODA(X0,Y0,tmin,tmax,nts,du_dt):
    t_span = (tmin,tmax)
    t = np.linspace(tmin,tmax,nts,endpoint = False)
    U0 = [X0,Y0]

    return solve_ivp(du_dt, t_span, U0,t_eval = t, method = 'LSODA')

#The volcotiy verlet (VelVerlet) integrator is configured to solve second order equations D_2(x) = f(t,x,v)
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

#Simple harmonic oscillator derivative (undamped). Take x=0 as the equillibrim point
def SHO(t,u): #we need to express the system using u. (we need t as an argument, but its not used)
    x,y = u #u is the phase space state of the system.

    #Put the derivative for x and y in here
    #x is the position. y is the velocity
    dx_dt = y
    dy_dt = -0.5*x
    return [dx_dt, dy_dt]

def SHO_damped(t,u): #we need to express the system using u. (we need t as an argument, but its not used)
    x,y = u #u is the phase space state of the system.

    #Put the derivative for x and y in here
    #x is the position. y is the velocity
    dx_dt = y
    dy_dt = -0.5*x  - 0.05*y
    return [dx_dt, dy_dt]

tmin = 0
tmax = 150
nts = 1500
sol = RK45(1,0,tmin,tmax,nts,SHO)

t = sol.t
X = sol.y[0]
Y = sol.y[1]

plt.plot(X,Y, label = "RK4(5)")

sol = LSODA(1,0,tmin,tmax,nts,SHO)

t = sol.t
X = sol.y[0]
Y = sol.y[1]

plt.plot(X,Y,label = "LSODA")

t, X,Y = VelVerlet(1,0,tmin,tmax,nts,SHO)
plt.plot(X,Y, label = "Velocity Verlet")
plt.legend()
plt.show()



sol = RK45(1,0,tmin,tmax,nts,SHO_damped)

t = sol.t
X = sol.y[0]
Y = sol.y[1]

plt.plot(X,Y, label = "RK4(5)")

sol = LSODA(1,0,tmin,tmax,nts,SHO_damped)

t = sol.t
X = sol.y[0]
Y = sol.y[1]

plt.plot(X,Y, label = "LSODA")

t, X,Y = VelVerlet(1,0,tmin,tmax,nts,SHO_damped)
plt.plot(X,Y, label = "Velocity Verlet")
plt.legend()
plt.show()