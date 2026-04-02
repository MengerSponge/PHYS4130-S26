import numpy as np
import matplotlib.pyplot as plt
from P3_Methods import Yoshida

# Simulation Parameters ----------------------------------------------------------------------
G = 1 #Gravitational Constant 

M =  [1.0, .50, 0.50, .450, 0.450]            # Masses
X0 = [[0,0], [2,0], [-2,0], [0,2], [0, -2]]          # Initial Positions
V0 = [[0,0], [0,1], [0, -1], [-1, 0], [1,0]]           # Initial Velocities

tmin = 0
tmax = 100
nts = 5000

# Funcitons to Compute Derivatives ----------------------------------------------------------------------
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

# Check That the Right Number of Conditions Have Been Given -----------------------------------------------------
if len(X0) != len(V0) or len(X0) != len(M):
    print()
    print("Inconsistent number of initial conditions or masses!")
    print("The following should all be equal:")
    print("Number of Point Masses:       ", len(M))
    print("Number of Initial Positions:  ", len(X0) )
    print("Number of Initial Velocities: ", len(V0))
    print()
    exit()

# Calculating and Plotting the Trajectories ------------------------------------------------------------------
    # Note: Strucutre of sol: 
    # sol = [times, [trajectory 1, velocity 1], [trajectory 2, velocity 2], ...]
    # where trajectories and velocities are two element np arrays.

sol = Yoshida(X0, V0, tmin, tmax, nts, Derivative) # All the parameters

for i in range(len(X0)):
    points = np.array(sol[i + 1][0])
    plt.plot(points[:, 0], points[:, 1]) #The : means to loop over all such elements

all_points = np.vstack([X[0] for X in sol[1:]])
plt.xlim(all_points[:, 0].min(), all_points[:, 0].max())
plt.ylim(all_points[:, 1].min(), all_points[:, 1].max())
plt.xlim(-7.5,7.5)
plt.ylim(-6,6)
plt.title("Trajectories")
plt.xlabel("X Position (m)")
plt.ylabel("Y Position (m)")
plt.savefig("5_Body_Trajectories.png", dpi = 400)
plt.show()

# Virial Theorem ------------------------------------------------------------------------------------------------------------------------------------
def Kinetic_Over_Time(sol): #takes the solution made by Yoshida method and finds the kinetic energy over time

    velocities = [V[1] for V in sol[1:]]
    T = [] #kinetic energy

    for i in range(len(sol[0])): #loop over the velocities at the current time
        v = [U[i] for U in velocities] #selects the ith term in velocities
        q = [0.5*M[j]*(np.linalg.norm(v[j])**2) for j in range(len(v))]

        T.append(sum(q))

    return T

def Potential(positions, masses, G): #takes in the positions at a time 

    positions = np.asarray(positions)
    masses = np.asarray(masses)

    N = len(masses)
    U = 0.0

    for i in range(N):
        for j in range(i + 1, N):
            r = np.linalg.norm(positions[i] - positions[j])
            if r != 0:
                U -= G * masses[i] * masses[j] / r

    return U

def Potential_Over_Time(sol):
    All_positions = [X[0] for X in sol[1:]]

    Potential_Energies = []
    for i in range(len(sol[0])): #do over time
        positions = [X[i] for X in All_positions] #goes over all the trajectoreis
        Potential_Energies.append(Potential(positions, M, G))

    return Potential_Energies


sol = Yoshida(X0, V0, tmin, tmax, nts, Derivative)

t = sol[0]
T = Kinetic_Over_Time(sol)
U = Potential_Over_Time(sol)

T_average = sum(T)/len(T)
U_average = sum(U)/len(U)
print(T_average)
print(U_average)

plt.plot(t, T, label = "Kinetic Energy")
plt.plot(t, U, label = "Potential Energy")
plt.plot(t, np.asarray(T) + np.asarray(U), label = "Total Energy")
plt.legend()
plt.title("Energy vs Time")
plt.xlabel("Time (s)")
plt.ylabel("Energy (J)")
plt.savefig("5_Body_Energies.png", dpi = 400)
plt.show()
