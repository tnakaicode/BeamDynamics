import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import time
from mpl_toolkits.mplot3d import axes3d

sys.path.append(os.path.join('../'))
from lib.ESPIC3D.esSolve import laplace2D
from lib.ESPIC3D.dirichlet import dirichlet as dirBC

NX = 50
NY = 75
LX = 1.25
LY = 2.3
DX = LX / NX
DY = LY / NY
X0 = 1.0
Y0 = 2.0

amplitude = 2.0


def nonGroundedWall(Yindex):
    return amplitude * np.sin(np.pi * Yindex / NY)


# Boundary conditions
V0x = dirBC(np.zeros((NY + 1)))
VNx = dirBC(np.fromfunction(nonGroundedWall, (NY + 1,)))

V0y = dirBC(np.zeros((NX + 1)))
VNy = dirBC(np.zeros((NX + 1)))

start = time.time()

solType = "direct"
solType = "jacobi"
solType ="gaussSeidel"

potential_1 = laplace2D(NX, DX, V0x, VNx, NY, DY, V0y, VNy,
                        solType, useCython=False, relTol=0.0, absTol=1.0)
potential_2 = laplace2D(NX, DX, V0x, VNx, NY, DY, V0y, VNy,
                        solType, useCython=False, relTol=0.0, absTol=0.5)
potential_3 = laplace2D(NX, DX, V0x, VNx, NY, DY, V0y, VNy,
                        solType, useCython=False, relTol=0.0, absTol=1.0e-3)

end = time.time()

print("That took", round(end - start, 1), "seconds.")


def plot2Darrays(arrays):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    if len(arrays) == 1:
        X = np.linspace(X0, X0 + LX, num=arrays[0].shape[0])
        Y = np.linspace(Y0, Y0 + LY, num=arrays[0].shape[1])
        X, Y = np.meshgrid(X, Y, indexing='ij')
    else:
        X = np.linspace(X0, X0 + LX, num=arrays[0][0].shape[0])
        Y = np.linspace(Y0, Y0 + LY, num=arrays[0][0].shape[1])
        X, Y = np.meshgrid(X, Y, indexing='ij')

    if len(arrays) == 1:
        ax.plot_wireframe(X, Y, arrays[0])
    else:
        for array in arrays:
            ax.plot_wireframe(
                X, Y, array[0][:], label='absolute tolerance = ' + array[1], color=array[2])

    ax.legend()


plot2Darrays([[potential_1, '1.0', 'green'],
              [potential_2, '0.5', 'red'],
              [potential_3, '1.0e-3', 'blue']])
plt.savefig("./ex_2d.png")

plot2Darrays([potential_1])
plt.savefig("./ex_2d_001.png")

plot2Darrays([potential_2])
plt.savefig("./ex_2d_002.png")

plot2Darrays([potential_3])
plt.savefig("./ex_2d_003.png")
