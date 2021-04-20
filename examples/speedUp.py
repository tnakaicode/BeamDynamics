import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import time
from mpl_toolkits.mplot3d import axes3d
sys.path.append(os.path.join('../'))

from solvers.dirichlet import dirichlet as dirBC
from solvers.esSolve import laplace2D

NX = 35
NY = 40
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
potential_NoCython = laplace2D(
    NX, DX, V0x, VNx, NY, DY, V0y, VNy, "gaussSeidel", relTol=0.0, absTol=1.0e-3, useCython=False)
end = time.time()
seconds_NoCython = end - start
print("No Cython took", round(seconds_NoCython, 1), "seconds.")

start = time.time()
potential_Cython = laplace2D(
    NX, DX, V0x, VNx, NY, DY, V0y, VNy, "gaussSeidel", relTol=0.0, absTol=1.0e-3, useCython=False)
end = time.time()
seconds_Cython = end - start
print("Cython took", round(seconds_Cython, 1), "seconds.")

print("Cython was", round(
    100.0 * (1.0 - seconds_Cython / seconds_NoCython), 1), "percent faster.")
