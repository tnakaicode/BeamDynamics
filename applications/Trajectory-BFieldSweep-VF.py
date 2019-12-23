import sys
sys.path.append('../lib/')
from BeamDynamicsTools import *
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
Rb = loadtxt('../data/CmodCoordinatesRZ.dat', usecols=[0])
Zb = loadtxt('../data/CmodCoordinatesRZ.dat', usecols=[1])

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
pl.figure(figsize=(20, 8))
for i in range(len(TrajectoryList)):
    pl.subplot(1, 2, 1)
    TrajectoryList[i].Plot2D('poloidal')
pl.subplot(1, 2, 1)
Vessel.Border('poloidal')
pl.xlim(0.2, 1.4)
pl.ylim(-0.7, 0.5)
pl.xlabel('R [m]')
pl.ylabel('Z [m]')
pl.title(r'Poloidal Projection ($\alpha$ = %0.1f$^o$, $\beta$ = %0.1f$^o$)' %
         (alpha0, beta0))
# pl.legend(Leg,loc=4)

# ------------------------------------------------------------------------------
# Plot 2D projections of Trajectories (Top View)
for i in range(len(TrajectoryList)):
    pl.subplot(1, 2, 2)
    TrajectoryList[i].Plot2D('top')
pl.subplot(1, 2, 2)
Vessel.Border('top')
pl.xlim(0, 1.2)
pl.ylim(-0.6, 0.6)
pl.xlabel('x [m]')
pl.ylabel('y [m]')
pl.title(r'Midplane Projection ($\alpha$ = %0.1f$^o$, $\beta$ = %0.1f$^o$)' %
         (alpha0, beta0))
ax = pl.subplot(1, 2, 2)
ax.legend(Leg, bbox_to_anchor=(1.28, 1.0))

# pl.legend(('B = 0.05','B = 0.10','B = 0.15','B = 0.20','B = 0.25','B = 0.30')

# ------------------------------------------------------------------------------
# Save Angular and Detection Quantities
if False:
    savetxt(OutputPath + 'geometry/TargetAngle_Vert_Horiz.dat', AngleComponents)
    savetxt(OutputPath + 'geometry/TargetCoordinates.dat', Coordinates)
    Header0 = '(0) I0 [A], (1) B0 [T], (2) X [m] , (3) Y [m], (4) Z [m], (5) incident angle [rad], (6) Detection Angle [rad], (7) optical path length [m] , (8) Detection Angle [rad], (9) Detection Angle [deg], (10) Detector Eff'
    savetxt(OutputPath + 'geometry/DetectionParameters.dat',
            (np.array(Parameters)), header=Header0)
# ------------------------------------------------------------------------------
# Save Figures
if True:
    FigName = 'TrajectoryProjections_alpha%2.2f_beta%2.2f' % (alpha0, beta0)
    FigPath = '/home/hbar/Dropbox/Research/AIMS/Magnet supply upgrade/Beam Modeling Results - Vertical Field Sweep/'
    pl.savefig(FigPath + FigName + '.pdf')
    pl.savefig(FigPath + FigName + '.png')
    print('File saved: ' + FigName)

pl.show()
