import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def Laplacian_9pt(U, h):
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

n = 3

U = np.array([[1,2,3],[4,5,6],[7,8,9]])

print(U)
print(np.roll(U, 1, 0))
print()
print(U)
print(np.roll(U, -1, 1))

U = np.array([[0, 0, 0],[0, 1, 0],[0, 0, 0]])
print(U)
print(Laplacian_9pt(U, 0.5))