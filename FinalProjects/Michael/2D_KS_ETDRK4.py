import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter

# Grid
N = 128
L = 32.0
t_steps = 1000
x = np.linspace(0, L, N, endpoint=False)
y = np.linspace(0, L, N, endpoint=False)
X, Y = np.meshgrid(x, y, indexing='ij')

# Fourier modes
kx = 2*np.pi * np.fft.fftfreq(N, d=L/N)
ky = 2*np.pi * np.fft.fftfreq(N, d=L/N)
Kx, Ky = np.meshgrid(kx, ky, indexing='ij')

K2 = Kx**2 + Ky**2

L = K2 - K2**2

# initial condtion and fourier transform of u0
u0 = 0.1*np.random.randn(N, N)
v0 = np.fft.fft2(u0)

# ETDRK4 coefficients ===================================
h = 0.25 #time step size

# operators 
E = np.exp(h*L)
E2 = np.exp(0.5*h*L)

# Phi functions compute with the contour trick by kassam
M = 16
r = np.exp(1j*np.pi*(np.arange(1, M+1) - 0.5)/M) #pieces of the unit circle

z = h * (L) #arguments of the phi functions
LR = z[:, :, None] + r[None, None, :] #take all the points in z and add a unit circle around them

# uses trapezoidal rule over the unit circle to approximate phi
phi1 = np.mean((np.exp(LR) - 1) / LR, axis=2)
phi2 = np.mean((np.exp(LR) - 1 - LR) / (LR**2), axis=2)
phi3 = np.mean((np.exp(LR) - 1 - LR - LR**2/2) / (LR**3), axis=2)

#phi1 = (np.exp(L) - 1) / L
#phi2 = (np.exp(L) - 1 - L) / (L**2)
#phi3 = (np.exp(L) - 1 - L - L**2/2) / (L**3)


print("Everything has been precomputed. Moving on to time-stepping.")
def NonLinear(u_ft): #takes in fourier transformed u and computes the non linear part.
    u_x = np.fft.ifft2(1j*Kx*u_ft).real
    u_y = np.fft.ifft2(1j*Ky*u_ft).real
    return np.fft.fft2(-0.5*(u_x**2 + u_y**2))

def Step(u_ft):
    # Stage a
    a = NonLinear(u_ft)

    # Stage b
    ua = E2 * (u_ft + 0.5*h*a)
    b = NonLinear(ua)

    # Stage c
    ub = E2 * (u_ft + 0.5*h*b)
    c = NonLinear(ub)

    # Stage d
    uc = E * (u_ft + h*c)
    d = NonLinear(uc)

    # Final update
    u_ft_next = E * u_ft + h * E * (phi1 * a + 2*phi2 * (b + c) + phi3 * d)
    return u_ft_next

v = v0
u_list = [u0]
t_list = [0]
nplt = 1 #save every nplt steps

for n in range(1,t_steps+1):
    t = n*h

    v = Step(v)
    if n%nplt == 0:
        u = np.fft.ifft2(v).real
        u_list.append(np.copy(u))
        t_list.append(np.copy(t))

fig, ax = plt.subplots()
im = ax.imshow(u_list[0], cmap = "plasma", interpolation='bilinear', extent=[0, h*(N-1), 0, h*(N-1)], origin='lower') #interpolation='bilinear'

fig.colorbar(im, ax = ax)

def update(frame):
    im.set_array(u_list[frame])
    return [im]

ani = FuncAnimation(fig, update, frames = len(t_list), interval = 1, blit=True, repeat = True)
print("Playing Animation.")
plt.show()