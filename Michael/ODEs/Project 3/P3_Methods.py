import numpy as np
from scipy.integrate import solve_ivp

#SciPy integrators --------------------------------------------------------------------------------
# These take in X0 and Y0 as lists of initial conditions
def RK45(X0, Y0, tmin, tmax, nts, du_dt):
    t_span = (tmin,tmax)
    t = np.linspace(tmin,tmax,nts,endpoint = True)
    
    X0 = np.asarray(X0)
    Y0 = np.asarray(Y0)

    solutions = [t]
    for i in range(len(X0)):
        sol = solve_ivp(du_dt, t_span, [X0[i], Y0[i]], t_eval = t, method = 'RK45')
        solutions.append([sol.y[0], sol.y[1]])

    return solutions

def LSODA(X0, Y0, tmin, tmax, nts, du_dt):
    t_span = (tmin,tmax)
    t = np.linspace(tmin,tmax,nts,endpoint = True)
    
    X0 = np.asarray(X0)
    Y0 = np.asarray(Y0)

    solutions = [t]
    for i in range(len(X0)):
        sol = solve_ivp(du_dt, t_span, [X0[i], Y0[i]], t_eval = t, method = 'LSODA')
        solutions.append([sol.y[0], sol.y[1]])

    return solutions

#Symplectic Integrators --------------------------------------------------------------------------------
# these take in X0 and Y0 as np arrays
# NOTE: These are written in a highly generalized form. They are configured to solve:
#
#           dx_dt = F(x,y,t) , dy_dt = G(x,y,t), -> du_dt = [F(x,y,t), G(x,y,t)]
#
# Which are not generally energy conserving.
# If you want to get energy conservation, you must express your system as:
#
#           dx_dt = y , dy_dt = a(x), -> du_dt = [y, a(x)]
#
# Where a is the acceleration.

# Velocity verlet (2nd order algorithm)
def VelVerlet(X0, Y0, tmin, tmax, nts, du_dt):
    t = np.linspace(tmin,tmax,nts,endpoint = True)

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

# 4th Order Yoshida Method (4th order version)
def Yoshida(X0, Y0, tmin, tmax, nts, du_dt):
    t = np.linspace(tmin,tmax,nts,endpoint = True)

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

#Analytic Solution --------------------------------------------------------------------------------
#Handles the different cases for under damped, over damped, and critically damped
def Analytic(X0, Y0, c, k, m, t):

    # Parameters:
    # X0 : initial position
    # Y0 : initial velocity
    # c  : damping coefficient
    # k  : spring constant
    # t  : time (can be scalar or numpy array)
    # m  : mass 

    # Natural frequency
    omega0 = np.sqrt(k / m)

    # Damping ratio
    gamma = c / (2 * m)

    # Discriminant
    disc = gamma**2 - omega0**2

    # --- Underdamped case ---
    if disc < 0:
        omega_d = np.sqrt(omega0**2 - gamma**2)

        A = X0
        B = (Y0 + gamma * X0) / omega_d

        return np.exp(-gamma * t) * (A * np.cos(omega_d * t) + B * np.sin(omega_d * t))

    # --- Critically damped case ---
    elif np.isclose(disc, 0): #Should reasonably handle issues with floating point calculation and comparison
        A = X0
        B = Y0 + gamma * X0

        return (A + B * t) * np.exp(-gamma * t)

    # --- Overdamped case ---
    else:
        r1 = -gamma + np.sqrt(disc)
        r2 = -gamma - np.sqrt(disc)

        # Solve for coefficients
        A = (Y0 - r2 * X0) / (r1 - r2)
        B = X0 - A

        return A * np.exp(r1 * t) + B * np.exp(r2 * t)

#Extra Methods I'll look at if I feel like it (probably not) --------------------------------------------------------------------------------
def Euler(X0, Y0, tmin, tmax, nts, f, g):
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

def RK2(X0, Y0, tmin, tmax, nts, f, g):


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
