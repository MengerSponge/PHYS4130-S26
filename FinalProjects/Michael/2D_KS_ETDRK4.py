import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter

# SIMULATION PARAMETERS ==================================================================================================================
#Control animation saving
save_anim = False

#color option for animation
color = 'viridis'

#nplt = number of steps stored in the animation. nplt = 1 means save all. larger means save less
nplt = 1

# N = number of points along a spatial dimension
N = 300
# size = length along the spatial dimensions
size = 30

#Tmax = maximum time
Tmax = 600

# h = time step size
h = 0.125 

#Verbosity controls how often the program reports progress
verbosity = int(Tmax/(10*h))

# initial condtions
x = size*np.arange(0,N)/N
y = size*np.arange(0,N)/N

# random noise ICS
u0 = 0.01*np.random.randn(len(x), len(y))

# single bump IC
u0 = np.zeros((len(x), len(y)))
u0[len(x)//2, len(y)//2] = 2

# CREATE THE SPATIAL GRID ==================================================================================================================

X, Y = np.meshgrid(x, y, indexing='ij')

# CREATE THE FOURIER MODES AND L OPERATOR ==================================================================================================================
dx = x[1] - x[0]
dy = y[1] - y[0]

kx = np.fft.fftfreq(N, d=dx) * 2*np.pi
ky = np.fft.fftfreq(N, d=dy) * 2*np.pi

Kx, Ky = np.meshgrid(kx, ky, indexing='ij')

K2 = Kx**2 + Ky**2

L = 0.5*(K2 - K2**2)

# ETDRK4 COEFFICIENTS =========================================================================================================
t_steps = int(Tmax/h)
# Operator Exponentials 
E = np.exp(h*L)
E2 = np.exp(0.5*h*L)

# Set up to do contour integral trick by Kasam
M = 20
r = np.exp(1j*np.pi*(np.arange(1, M+1) - 0.5)/M) #pieces of the unit circle
z = h * (L) #arguments of the phi functions
LR = z[:, :, None] + r[None, None, :] #take all the points in z and add a unit circle around them

# the actual coefficents computed using the contour
f1 = h*np.mean(((-4-LR+np.exp(LR)*(4-3*LR+LR**2))/LR**3),axis=2).real
f2 = h*np.mean(((2+LR+np.exp(LR)*(-2+LR))/LR**3),axis=2).real
f3 = h*np.mean(((-4-3*LR-LR**2+np.exp(LR)*(4-LR))/LR**3),axis=2).real

Q  = h*np.mean(((np.exp(LR/2) - 1) / (LR)), axis=2).real

print("Everything has been precomputed. Moving on to time-stepping.")

# TIME STEPPING =========================================================================================================
def NonLinear(u_ft): #takes in fourier transformed u and computes the non linear part.
    #spatial derivatives
    u_x = np.fft.ifft2(-1j*Kx*u_ft).real
    u_y = np.fft.ifft2(-1j*Ky*u_ft).real

    # nonlinear terms in position space
    N = -0.5 * (u_x**2 + u_y**2)
    return np.fft.fft2(N) #return to fourier space

def Step(u_ft):
    a00 = u_ft[0,0] #to enforce constant average
    t = n*h

    #RK4 stepping
    Nu = NonLinear(u_ft)
    a = E2*u_ft + Q*Nu
    Na = NonLinear(a)
    b = E2*u_ft + Q*Na
    Nb = NonLinear(b)
    c = E2*a + Q*(2*Nb-Nu)
    Nc = NonLinear(c)
    u_ft_new = E*u_ft + Nu*f1 + 2*(Na+Nb)*f2 + Nc*f3

    #enforce constant average
    u_ft_new[0,0] = a00
    return u_ft_new

v = np.fft.fft2(u0) #create v as the fourier transform of u0

#start these lists to save the funciton over time
u_list = [u0]
t_list = [0]

for n in range(1,t_steps+1):
    #step
    t = n*h
    v = Step(v)

    if n%nplt == 0: #to save for animating
        u = np.fft.ifft2(v).real
        u_list.append(np.copy(u))
        t_list.append(np.copy(t))
    if n%verbosity == 0: # to output progress
        print(100*n/len(range(1,t_steps+1)), " % complete")
    
# ANIMATING =========================================================================================================
print("Building animation")
#max an min for setting plot range
umax = np.max(np.array(u_list))
umin = np.min(np.array(u_list))
print("Maximum: ", umax, " Minimum: ", umin)

# Animation junk
fig, ax = plt.subplots()
plt.title("2D Kuramoto-Sivashinksy Equation with ETDRK4", y = 1.05)
ax.set_xlabel("x")
ax.set_ylabel("y")
im = ax.imshow(u_list[0], cmap = color, interpolation='bilinear', extent=[0, size, 0, size], origin='lower', vmin = umin, vmax = umax) #interpolation='bilinear'

fig.colorbar(im, ax = ax)

def update(frame):
    im.set_array(u_list[frame])
    return [im]

ani = FuncAnimation(fig, update, frames = len(t_list), interval = 1, blit=True, repeat = True)

#different choice for the animation
if save_anim == True:
    print("Saving Animation.")
    writervideo = FFMpegWriter(fps=60)
    ani.save('PDE_Solution.mp4', writer=writervideo)
else:
    print("Playing Animation.")
    plt.show()
