import numpy as np
import matplotlib.pyplot as plt

from concave_hull import concave_hull
from shapely.geometry import Polygon
from P3_Methods import VelVerlet, Yoshida, LSODA

#constants ---------------------------------------
k = 1.2 # N/m        (spring constant)
m = 1 # kg         (mass)
c = 0.25 # N/(m/s) (damping term strength)

tmin = 0 #s start time
tmax = 150 #s end time
nts = 1500 #number of points between tmin and tmax

#derivatives
def SHO(t, u): #taken from P3_SHO.py
    x,y = u

    dx_dt = y
    dy_dt = -(k/m)*x
    return [dx_dt, dy_dt]

def SHO_damped(t, u): 
    x,y = u

    #Put the derivative for x and y in here
    #x is the position. y is the velocity
    dx_dt = y
    dy_dt = -(k/m)*x  - (c/m)*y
    return [dx_dt, dy_dt]

#Arrays of the initial conditions. These should work with my methods
#since numpy attributes are vectorized
n = 4
X0 = np.random.uniform(0.99, 1, size = n)
Y0 = np.random.uniform(0.99, 1, size = n)

solutions = VelVerlet(X0, Y0, tmin, tmax, nts, SHO)
t = solutions[0]

for i in range(1,len(X0)):
    plt.plot(solutions[i][0], solutions[i][1], label = "Initial Conditions "+str(i))

plt.title("Phase Space")
plt.xlabel("Position (m)")
plt.ylabel("Velocity (m/s)")
plt.show()

PhaseVolumes = np.zeros(len(t))
for j in range(len(t)):
    P = [[solutions[k][0][j], solutions[k][1][j]] for k in range(1,n)]
    P = concave_hull(P, concavity=2.0)
    poly = Polygon(P)
    PhaseVolumes[j] = poly.area

plt.plot(t, PhaseVolumes, label = "Velocity Verlet")

solutions = Yoshida(X0, Y0, tmin, tmax, nts, SHO)
t = solutions[0]

PhaseVolumes = np.zeros(len(t))
for j in range(len(t)):
    P = [[solutions[k][0][j], solutions[k][1][j]] for k in range(1,n)]
    P = concave_hull(P, concavity=2.0)
    poly = Polygon(P)
    PhaseVolumes[j] = poly.area

plt.plot(t, PhaseVolumes, label = "4th Order Yoshida")

solutions= LSODA(X0, Y0, tmin, tmax, nts, SHO)
t = solutions[0]

PhaseVolumes = np.zeros(len(t))
for j in range(len(t)):
    P = [[solutions[k][0][j], solutions[k][1][j]] for k in range(1,n)]
    P = concave_hull(P, concavity=2.0)
    poly = Polygon(P)
    PhaseVolumes[j] = poly.area

plt.plot(t, PhaseVolumes, label = "LSODA")
plt.title("Phase Space Volume vs Time")
plt.ylabel("Phase Space Volume (m * m/s)")
plt.xlabel("Time (s)")
plt.legend()
plt.ylim(0, 5e-5)
plt.savefig("SHO_Volume.png", dpi = 400)
plt.show()




solutions = VelVerlet(X0, Y0, tmin, tmax, nts, SHO_damped)
t = solutions[0]

for i in range(1,len(X0)):
    plt.plot(solutions[i][0], solutions[i][1], label = "Initial Conditions "+str(i))

plt.title("Phase Space")
plt.xlabel("Position (m)")
plt.ylabel("Velocity (m/s)")
plt.show()

PhaseVolumes = np.zeros(len(t))
for j in range(len(t)):
    P = [[solutions[k][0][j], solutions[k][1][j]] for k in range(1,n)]
    P = concave_hull(P, concavity=2.0)
    poly = Polygon(P)
    PhaseVolumes[j] = poly.area

plt.plot(t, PhaseVolumes, label = "Velocity Verlet")

solutions = Yoshida(X0, Y0, tmin, tmax, nts, SHO_damped)
t = solutions[0]

PhaseVolumes = np.zeros(len(t))
for j in range(len(t)):
    P = [[solutions[k][0][j], solutions[k][1][j]] for k in range(1,n)]
    P = concave_hull(P, concavity=2.0)
    poly = Polygon(P)
    PhaseVolumes[j] = poly.area

plt.plot(t, PhaseVolumes, label = "4th Order Yoshida")

solutions= LSODA(X0, Y0, tmin, tmax, nts, SHO_damped)
t = solutions[0]

PhaseVolumes = np.zeros(len(t))
for j in range(len(t)):
    P = [[solutions[k][0][j], solutions[k][1][j]] for k in range(1,n)]
    P = concave_hull(P, concavity=2.0)
    poly = Polygon(P)
    PhaseVolumes[j] = poly.area

plt.plot(t, PhaseVolumes, label = "LSODA")
plt.title("Phase Space Volume vs Time")
plt.ylabel("Phase Space Volume (m * m/s)")
plt.xlabel("Time (s)")
plt.legend()
#plt.ylim(0, 0.0005)
plt.savefig("SHO_Damped_Volume.png", dpi = 400)
plt.show()