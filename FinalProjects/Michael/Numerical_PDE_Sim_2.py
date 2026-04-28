import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.signal import convolve2d
from matplotlib.animation import FuncAnimation, FFMpegWriter
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

# Functions for creating Derivatives and Different Fixed Points ---------------------------------------------------------
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

def Partial_x_1(U, h):
    return(-np.roll(U, -2, 1) + 8*np.roll(U, -1, 1) - 8*np.roll(U, 1, 1) + np.roll(U, 2, 1))/(12*h)

def Partial_y_1(U, h):
    return (-np.roll(U, -2, 0) + 8*np.roll(U, -1, 0) - 8*np.roll(U, 1, 0) + np.roll(U, 2, 0))/(12*h)

def BlockBoundary(x0, x1, y0, y1, N): #produces the mask for a fixed block points
    mask = np.zeros((N,N), dtype = bool)
    mask[y0:y1, x0:x1] = 1
    return mask

# System ---------------------------------------------------------------------
Save_Animation = False

N = 151 # Num points on a side
h = 0.5 # spatial step size

T_final = 100
Nt = 400

dt = T_final/Nt

# Initial conditions
U0 = np.zeros((N, N))
#U0[13:15, N//2 - 50:N//2-30] = 0
U0[120:130, N//2 - 5:N//2 ] = 20
V0 = np.zeros((N, N))

# Package the full system as one long array
Y0 = np.concatenate([U0.flatten(), V0.flatten()])

# Fixed Point Locations (values must be altered in the rhs function) 1 = fix this point
Do_BCS = True
left_boundary   = BlockBoundary(0  , 1, 0  , N, N)
right_boundary  = BlockBoundary(N-1, N, 0  , N, N)
top_boundary    = BlockBoundary(0  , N, 0  , 1, N)
bottom_boundary = BlockBoundary(0  , N, N-1, N, N)

slit = BlockBoundary(0, N//2 - 6, N//2-20, N//2 -14, N) + BlockBoundary(N//2 + 6, N, N//2-20, N//2 -14, N) + BlockBoundary(N//2 - 1, N//2 + 1,N//2-20, N//2 -14, N)
#System Derivatives
alpha = 2
def F(t, U, V):
    dU_dt = V
    dV_dt = alpha**2 *Laplacian_9pt(U, h) -0.01*V

    return [dU_dt, dV_dt]

# Package up the system for ODEint
def rhs(t, Y):
    
    U = Y[:N*N].reshape((N, N)) #slice up to N*N - 1
    V = Y[N*N:].reshape((N, N)) #slice ater N*N
    
    #Enforce values
    if Do_BCS == True:
        U[left_boundary]   = 0
        U[right_boundary ] = 0
        U[top_boundary]    = 0
        U[bottom_boundary] = 0
        U[slit] = 0

        V[left_boundary]   = 0
        V[right_boundary ] = 0 
        V[top_boundary]    = 0
        V[bottom_boundary] = 0
        V[slit] = 0
    
    # compute derivatives...
    dU_dt, dV_dt = F(t, U, V)

    #Enforce dynamics
    if Do_BCS == True:
        dU_dt[left_boundary]   = 0
        dU_dt[right_boundary ] = 0
        dU_dt[top_boundary]    = 0
        dU_dt[bottom_boundary] = 0
        dU_dt[slit] = 0

        dV_dt[left_boundary]   = 0
        dV_dt[right_boundary ] = 0 
        dV_dt[top_boundary]    = 0
        dV_dt[bottom_boundary] = 0
        dV_dt[slit] = 0

    return np.concatenate([dU_dt.flatten(), dV_dt.flatten()])

# Solving the system ----------------------------------------------------------
print("Computing Solution.")
sol = solve_ivp(
    rhs, 
    t_span = (0, T_final), 
    y0 = Y0,
    method = 'RK45',
    t_eval = np.linspace(0, T_final, Nt)
    )
# solution over time

Uframes = [5*sol.y[:N*N,k].reshape((N,N)) for k in range(0, Nt)]
Vframes = [sol.y[N*N:,k].reshape((N,N)) for k in range(0, Nt)]

print("Solution Computed. Building Animation.")

# Animation ----------------------------------------------------------------------

custom_cmap = LinearSegmentedColormap.from_list(
    "blue_black_red",
    [(0, "blue"), (0.5, "black"), (1, "red")])

custom_norm = TwoSlopeNorm(
    vmin=-5,
    vcenter=0,
    vmax=5)

fig, ax = plt.subplots()
im = ax.imshow(Uframes[0], cmap = custom_cmap, norm = custom_norm, interpolation='bilinear', extent=[0, h*(N-1), 0, h*(N-1)], origin='lower') #interpolation='bilinear'
#title = ax.set_title(f"t = {sol.t[0]:.3f}")

fig.colorbar(im, ax = ax)

def update(frame):
    im.set_array(Uframes[frame])
    #title.set_text(f"t = {sol.t[frame]:.3f}")
    return [im]

ani = FuncAnimation(fig, update, frames = Nt, interval = 1, blit=True, repeat = True)
print("Playing Animation.")
plt.show()

if Save_Animation == True:
    print("Saving Animation")
    ani = FuncAnimation(fig, update, frames = Nt, interval = 1, blit=True, repeat = True)
    writervideo = FFMpegWriter(fps=60)
    ani.save('PDE_Solution.mp4', writer=writervideo)
    plt.close()
    print("Animation Saved as PDE_Solution.mp4")