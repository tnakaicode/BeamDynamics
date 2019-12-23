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

import numpy as np
import pylab as pl

S1 = np.matrix([
    [0.577100, 0.398000, 0.000000, 0.000000, 0.000000, 0.000000],
    [0.398000, 171.8262, 0.000000, 0.000000, 0.000000, 0.000000],
    [0.000000, 0.000000, 0.343900, -0.27150, 0.000000, 0.000000],
    [0.000000, 0.000000, -0.27150, 238.3722, 0.000000, 0.000000],
    [0.000000, 0.000000, 0.000000, 0.000000, 1.297156, 2.343722],
    [0.000000, 0.000000, 0.000000, 0.000000, 2.343722, 134.9344]], float)

S0 = np.matrix([
    [1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000],
    [0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000],
    [0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000],
    [0.000000, 0.000000, 0.000000, 1.000000, 0.000000, 0.000000],
    [0.000000, 0.000000, 0.000000, 0.000000, 1.000000, 0.000000],
    [0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1.000000]], float)


E0 = Ellipse(S0)

E0.PlotProjectionXY()
E0.PlotProjectionXY(Axz=np.pi / 4, Ayz=0 * np.pi / 3)


w = 2
plt.xlim(-w, w)
plt.ylim(-w, w)
plt.show()
