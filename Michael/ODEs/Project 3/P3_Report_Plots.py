import numpy as np
import matplotlib.pyplot as plt
from P3_Methods import RK45, LSODA, VelVerlet

def SHO(t, u): #we need to express the system using u = [x, y]. (we need t as an argument, but its not always used)
    x,y = u #u is the phase space state of the system.

    #Put the derivative for x and y in here
    #x is the position. y is the velocity
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

def Energy(X,Y): #Energy at different values of time
    return  0.5*m*(Y**2) + 0.5*k*(X**2)


k = 1 # N/m        (spring constant)
m = 1 # kg         (mass)
c = 0.5 # N/(m/s) (damping term strength)

tmin = 0 #s start time
tmax = 50 #s end time
nts = 100 #number of points between tmin and tmax

X0 = 1 #m Initial positon
Y0 = 0 #m/s Initial velocity

E0 = 0.5*m*(Y0**2) + 0.5*k*(X0**2)


# FIGURE 1
Energies = []

solution = RK45([X0], [Y0], tmin, tmax, nts, SHO)
XRK_Undamped = solution[1][0]
YRK_Undamped = solution[1][1]
Energies.append(Energy(XRK_Undamped, YRK_Undamped))

plt.plot(XRK_Undamped,YRK_Undamped, label = "RK4(5)")


solution = LSODA([X0], [Y0], tmin, tmax, nts, SHO)
XLSODA_Undamped = solution[1][0] 
YLSODA_Undamped = solution[1][1]
Energies.append(Energy(XLSODA_Undamped, YLSODA_Undamped))

plt.plot(XLSODA_Undamped,YLSODA_Undamped, label = "LSODA")


solution = VelVerlet([X0], [Y0], tmin, tmax, nts, SHO)
XVerlet_Undamped = solution[1][0]
YVerlet_Undamped = solution[1][1]
Energies.append(Energy(XVerlet_Undamped, YVerlet_Undamped))

plt.plot(XVerlet_Undamped, YVerlet_Undamped, label = "Velocity Verlet")


plt.title("Phase Space Plots (Undamped)")
plt.xlabel("Position (m)")
plt.ylabel("Velocity (m/s)")
plt.legend(loc = 'upper right')
plt.savefig("SHO_Phase.png",dpi = 400)
plt.show()

#FIGURE 2
fig, ax = plt.subplots(1, 3, figsize=(15 ,5))
t = solution[0]
ax[0].plot(t, Energies[0], label = "RK4(5)", color = "blue")
ax[0].axhline(y = E0, color = 'r', linestyle = '--', label = "True Energy")
ax[0].legend(loc = 'upper right')
ax[0].set_xlabel("Time (s)")
ax[0].set_ylabel("Energy (J)")



ax[1].plot(t, Energies[1], label = "LSODA", color = "orange")
ax[1].axhline(y = E0, color = 'r', linestyle = '--', label = "True Energy")
ax[1].legend(loc = 'upper right')
ax[1].set_xlabel("Time (s)")
ax[1].set_ylabel("Energy (J)")

ax[2].plot(t, Energies[2], label = "Velocity Verlet", color = "green")
ax[2].axhline(y = E0, color = 'r', linestyle = '--', label = "True Energy")
ax[2].legend(loc = 'upper right')
ax[2].set_xlabel("Time (s)")
ax[2].set_ylabel("Energy (J)")

plt.tight_layout()
plt.savefig("SHO_Energy.png",dpi = 400)
plt.show()

# Figure 3

k = 1 # N/m        (spring constant)
m = 1 # kg         (mass)
c = 0.5 # N/(m/s) (damping term strength)

tmin = 0 #s start time
tmax = 50 #s end time
nts = 100 #number of points between tmin and tmax

X0 = 1 #m Initial positon
Y0 = 0 #m/s Initial velocity

E0 = 0.5*m*(Y0**2) + 0.5*k*(X0**2)

Energies = []

solution = RK45([X0], [Y0], tmin, tmax, nts, SHO_damped)
XRK_Undamped = solution[1][0]
YRK_Undamped = solution[1][1]
Energies.append(Energy(XRK_Undamped, YRK_Undamped))

plt.plot(XRK_Undamped,YRK_Undamped, label = "RK4(5)")


solution = LSODA([X0], [Y0], tmin, tmax, nts, SHO_damped)
XLSODA_Undamped = solution[1][0] 
YLSODA_Undamped = solution[1][1]
Energies.append(Energy(XLSODA_Undamped, YLSODA_Undamped))

plt.plot(XLSODA_Undamped,YLSODA_Undamped, label = "LSODA")


solution = VelVerlet([X0], [Y0], tmin, tmax, nts, SHO_damped)
XVerlet_Undamped = solution[1][0]
YVerlet_Undamped = solution[1][1]
Energies.append(Energy(XVerlet_Undamped, YVerlet_Undamped))

plt.plot(XVerlet_Undamped, YVerlet_Undamped, label = "Velocity Verlet")


plt.title("Phase Space Plots (Damped)")
plt.xlabel("Position (m)")
plt.ylabel("Velocity (m/s)")
plt.legend(loc = 'upper right')
plt.savefig("SHO_Damped_Phase.png",dpi = 400)
plt.show()

#FIGURE 4
fig, ax = plt.subplots(1, 3, figsize=(15 ,5))
t = solution[0]
ax[0].plot(t, Energies[0], label = "RK4(5)", color = "blue")
ax[0].axhline(y = E0, color = 'r', linestyle = '--', label = "Initial Energy")
ax[0].legend(loc = 'upper right')
ax[0].set_xlabel("Time (s)")
ax[0].set_ylabel("Energy (J)")



ax[1].plot(t, Energies[1], label = "LSODA", color = "orange")
ax[1].axhline(y = E0, color = 'r', linestyle = '--', label = "Initial Energy")
ax[1].legend(loc = 'upper right')
ax[1].set_xlabel("Time (s)")
ax[1].set_ylabel("Energy (J)")

ax[2].plot(t, Energies[2], label = "Velocity Verlet", color = "green")
ax[2].axhline(y = E0, color = 'r', linestyle = '--', label = "Initial Energy")
ax[2].legend(loc = 'upper right')
ax[2].set_xlabel("Time (s)")
ax[2].set_ylabel("Energy (J)")

plt.tight_layout()
plt.savefig("SHO_Damped_Energy.png",dpi = 400)
plt.show()