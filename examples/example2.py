import numpy as np
import matplotlib.pyplot as plt
import sys
import os
sys.path.append(os.path.join('../'))
from lib.BeamDynamicsTools.Boundary import Boundary
from lib.BeamDynamicsTools.Bfield import Bfield, BfieldTF, BfieldVF
from lib.BeamDynamicsTools.Trajectory import Trajectory
from lib.BeamDynamicsTools.Beam import Beam

# ===============================================================================
# Calculate trajectory matrix evolution for 4 values of toroidal B
# ===============================================================================

# ------------------------------------------------------------------------------
# Import poloidal boundary points
Rb = np.loadtxt('../data/CmodCoordinatesRZ.dat', usecols=[0])
Zb = np.loadtxt('../data/CmodCoordinatesRZ.dat', usecols=[1])

# ------------------------------------------------------------------------------
# Generate vessel Boundary
Vessel = Boundary(Rb, Zb)

# ------------------------------------------------------------------------------
# 3D plot of vessel boundary
ax = Vessel.Figure3D(1)
Vessel.Plot3D(ax)

# ------------------------------------------------------------------------------
# Inputs for four B-field settings
In = np.array([0.0, 1600.0, 3120, 4450.0])
Bn = np.array([0.0, 0.05818182, 0.11345455, 0.16181818])

# ===============================================================================
# Perform Trajectory and sigma dynamics calculation for B-Field Sweep
# ===============================================================================

AngleComponents = []
Coordinates = []
Parameters = []
TrajectoryList = []
OutputPath = './output/'
for i in [0, 1, 2, 3]:  # range(len(Bn)):
    B = BfieldTF(B0=Bn[i])
    Bv = BfieldVF(B0=0.00000)
    T = Trajectory(Vessel, B, Bv)
    TrajectoryList.append(T)

# ------------------------------------------------------------------------------
# Save target parameters
#	T.target.SaveTargetParameters(TFCurrent=In[i],Path=OutputPath+'geometry/')

# ------------------------------------------------------------------------------
# append lists of Target Quantities
    AngleComponents.append([T.target.VAngle, T.target.HAngle])
    Coordinates.append([T.target.R, T.target.Z, T.target.Phi])
    Parameters.append(T.target.GetDetectionParameters())

# ------------------------------------------------------------------------------
# Plot 3D trajectory results
Color = ['b', 'g', 'r', 'c']
for i in range(len(TrajectoryList)):
    TrajectoryList[i].LineColor = Color[i]
    TrajectoryList[i].Plot3D(ax)
#	TrajectoryList[i].target.Plot3D(ax);

TrajectoryList[-1].Limits3D(ax)

# Plot 2D projections of Trajectories
#	plt.figure(10); T.Plot2D()
#	plt.figure(11); T.Plot2D('top')
#plt.figure(10); Vessel.Border(); plt.xlim(0.2,1.4); plt.ylim(-0.7,0.5)
#plt.xlabel('R [m]'); plt.ylabel('Z [m]'); plt.title('Poloidal Projection')
#plt.figure(11); Vessel.Border('top'); plt.xlim(0,1.2); plt.ylim(-0.6,0.6)
#plt.xlabel('x [m]'); plt.ylabel('y [m]'); plt.title('Midplane Projection')


# Save Angular and Detection Quantities
if False:
    np.savetxt(OutputPath + 'geometry/TargetAngle_Vert_Horiz.dat',
               AngleComponents)
    np.savetxt(OutputPath + 'geometry/TargetCoordinates.dat', Coordinates)
    Header0 = '(0) I0 [A], (1) B0 [T], (2) X [m] , (3) Y [m], (4) Z [m], (5) incident angle [rad], (6) Detection Angle [rad], (7) optical path length [m] , (8) Detection Angle [rad], (9) Detection Angle [deg], (10) Detector Eff'
    np.savetxt(OutputPath + 'geometry/DetectionParameters.dat',
               (np.np.array(Parameters)), header=Header0)

plt.show()
