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
import pylab as pl
import matplotlib as mpl

# ===============================================================================
# Calculates beam trajectories over a parameter sweep over a range of toroidal
# field settings
# ===============================================================================

# ------------------------------------------------------------------------------
# Define np.array of injection angles
# (x,y,z) = (1.798m, -0.052m, 0.243m)
#  alpha = 12.6 degrees (X-Z plane)
#  beta = 8.0 degrees (X-Y plane)

alpha0 = 12.6
beta0 = 8.0

alpha = alpha0 / 180.0 * np.pi
beta = beta0 / 180.0 * np.pi
print(alpha, beta)
Rinjection = [1.798, -0.052, 0.243]
Vinjection = [-np.cos(alpha) * np.cos(beta), np.cos(alpha) * np.sin(beta), -np.sin(alpha)]
#Energy = [0.594e6, 0.740e6, 0.900e6]
Energy = 0.9e6  # np.linspace(0.594e6,0.900e6,10)

# ------------------------------------------------------------------------------
# Import poloidal Boundary points
Rb = np.loadtxt('../data/CmodCoordinatesRZ.dat', usecols=[0])
Zb = np.loadtxt('../data/CmodCoordinatesRZ.dat', usecols=[1])

# ------------------------------------------------------------------------------
# Generate vessel Boundary
Vessel = Boundary(Rb, Zb)

# ------------------------------------------------------------------------------
# 3D plot of vessel Boundary
ax = Vessel.Figure3D()
Vessel.Plot3D(ax)

# ------------------------------------------------------------------------------
# Inputs for four B-field settings
Bn = np.linspace(-0.05, 0.05, 21)

# ------------------------------------------------------------------------------
# Generate Color Map
CMAP = mpl.colors.LinearSegmentedColormap.from_list(
    'mycolors', ['green', 'blue', 'black', 'red', 'orange'])

# ===============================================================================
# Perform Trajectory calculation for B-Field Sweep
# ===============================================================================

AngleComponents = []
Coordinates = []
Parameters = []
TrajectoryList = []
OutputPath = '../output/'
# Color=['k','g','r','c','b','m','g','r','c','b','m','g']

for i in range(len(Bn)):
    B = BfieldTF(B0=0.0)
    Bv = BfieldVF(B0=Bn[i])
    T = Trajectory(Vessel, B, Bv, v0=Vinjection, T0=Energy)
    T.LineColor = CMAP(1.0 * i / len(Bn))
    T.LineWidth = 2.0
    TrajectoryList.append(T)

# ------------------------------------------------------------------------------
    # Save Target parameters
#	T.Target.SaveTargetParameters(TFCurrent=In[i],Path=OutputPath+'geometry/')

    # append lists of Target Quantities
#	AngleComponents.append([T.Target.VAngle,T.Target.HAngle])
#	Coordinates.append([T.Target.R,T.Target.Z,T.Target.Phi])
#	Parameters.append(T.Target.GetDetectionParameters())

# ------------------------------------------------------------------------------
# Plot 3D results

for i in range(len(TrajectoryList)):
    TrajectoryList[i].Plot3D(ax)
    TrajectoryList[i].target.Plot3D(ax)
TrajectoryList[-1].Limits3D(ax)

# ------------------------------------------------------------------------------
# Construct Legend
Leg = []
for i in range(len(Bn)):
    Leg.append('NI = %0.1fkA' % (TrajectoryList[i].BFieldVF.I0 / 1e3))

# ------------------------------------------------------------------------------
# Plot 2D projections of Trajectories (Poloidal View)
plt.figure(figsize=(20, 8))
for i in range(len(TrajectoryList)):
    plt.subplot(1, 2, 1)
    TrajectoryList[i].Plot2D('poloidal')
plt.subplot(1, 2, 1)
Vessel.Border('poloidal')
plt.xlim(0.2, 1.4)
plt.ylim(-0.7, 0.5)
plt.xlabel('R [m]')
plt.ylabel('Z [m]')
plt.title(r'Poloidal Projection ($\alpha$ = %0.1f$^o$, $\beta$ = %0.1f$^o$)' %
         (alpha0, beta0))
# plt.legend(Leg,loc=4)

# ------------------------------------------------------------------------------
# Plot 2D projections of Trajectories (Top View)
for i in range(len(TrajectoryList)):
    plt.subplot(1, 2, 2)
    TrajectoryList[i].Plot2D('top')
plt.subplot(1, 2, 2)
Vessel.Border('top')
plt.xlim(0, 1.2)
plt.ylim(-0.6, 0.6)
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.title(r'Midplane Projection ($\alpha$ = %0.1f$^o$, $\beta$ = %0.1f$^o$)' %
         (alpha0, beta0))
ax = plt.subplot(1, 2, 2)
ax.legend(Leg, bbox_to_anchor=(1.28, 1.0))

# plt.legend(('B = 0.05','B = 0.10','B = 0.15','B = 0.20','B = 0.25','B = 0.30')

# ------------------------------------------------------------------------------
# Save Angular and Detection Quantities
if False:
    np.savetxt(OutputPath + 'geometry/TargetAngle_Vert_Horiz.dat', AngleComponents)
    np.savetxt(OutputPath + 'geometry/TargetCoordinates.dat', Coordinates)
    Header0 = '(0) I0 [A], (1) B0 [T], (2) X [m] , (3) Y [m], (4) Z [m], (5) incident angle [rad], (6) Detection Angle [rad], (7) optical path length [m] , (8) Detection Angle [rad], (9) Detection Angle [deg], (10) Detector Eff'
    np.savetxt(OutputPath + 'geometry/DetectionParameters.dat',
            (np.array(Parameters)), header=Header0)
# ------------------------------------------------------------------------------
# Save Figures
if True:
    FigName = 'TrajectoryProjections_alpha%2.2f_beta%2.2f' % (alpha0, beta0)
    FigPath = '/home/hbar/Dropbox/Research/AIMS/Magnet supply upgrade/Beam Modeling Results - Vertical Field Sweep/'
    plt.savefig(FigPath + FigName + '.pdf')
    plt.savefig(FigPath + FigName + '.png')
    print('File saved: ' + FigName)

plt.show()
