import matplotlib.pyplot as plt
import os
import sys
sys.path.append(os.path.join('../'))
#sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '/../solvers')

from solvers.dirichlet import dirichlet as dirBC
from solvers.particle import particle
from solvers.particleUtils import velocityToMomentum, momentumToVelocity
from solvers.esSolve import laplace1D, potentialToElectricField, electricFieldAtPoint


NX = 100
LX = 1.2

DX = LX / NX

V0 = dirBC(1.0)
VN = dirBC(2.0)

X0 = 1.5

pot1D = laplace1D(NX, DX, V0, VN, "gaussSeidel", relTol=0.0,
                  absTol=1.0e-3, useCython=False)

# could use X0 here to plot versus position instead of index
plt.plot(pot1D)
plt.savefig("./ex_1d.png")
