import sys
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
sys.path.append(os.path.join('../'))
from lib.BeamDynamicsTools.Boundary import Boundary
from lib.BeamDynamicsTools.Bfield import Bfield, BfieldTF, BfieldVF
from lib.BeamDynamicsTools.Trajectory import Trajectory
from lib.BeamDynamicsTools.Beam import Beam
from lib.BeamDynamicsTools.Ellipse import Ellipse
from BoundaryStellarator import *
import numpy as np
import pylab as pl
import timeit

# Import poloidal boundary points
Rb = np.array(np.loadtxt('../data/CmodCoordinatesRZ.dat', usecols=[0]))  # -0.66
Zb = np.loadtxt('../data/CmodCoordinatesRZ.dat', usecols=[1])

# ------------------------------------------------------------------------------
# Generate vessel boundary
Vessel = BoundaryStellarator(Rb, Zb)

# ------------------------------------------------------------------------------
# 3D plot of vessel boundary
ax = Vessel.Figure3D(3)
Vessel.Plot3D(ax)


#
plt.figure(1)
Vessel.Border()

plt.figure(2)
Vessel.Plot2D(2)
