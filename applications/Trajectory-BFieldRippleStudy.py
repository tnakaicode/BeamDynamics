import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sp
import sys
import os
from numpy.linalg import norm, det
sys.path.append(os.path.join('../'))
from lib.BeamDynamicsTools.Boundary import Boundary
from lib.BeamDynamicsTools.Bfield import Bfield, BfieldTF, BfieldVF, CalculateB0, CalculateI0
from lib.BeamDynamicsTools.Trajectory import Trajectory
from lib.BeamDynamicsTools.Beam import Beam
from lib.BeamDynamicsTools.Ellipse import Ellipse
import pylab as pl
import matplotlib as mpl

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
Bn = np.linspace(-0.45, 0.45, 41)
#Bn = np.array([0.0])

# ------------------------------------------------------------------------------
# Generate Color Map
#CMAP = mpl.colors.LinearSegmentedColormap.from_list('mycolors',['black','red','orange'])
CMAP = mpl.colors.LinearSegmentedColormap.from_list(
    'mycolors', ['green', 'blue', 'black', 'red', 'orange'])

AngleComponents = []
Coordinates = []
Parameters = []
TrajectoryList = []
OutputPath = '../output/'
# Color=['k','g','r','c','b','m','g','r','c','b','m','g']

# ===============================================================================
# Calculate Target Error vs Magnet Ripple
# ===============================================================================

# ------------------------------------------------------------------------------
# Define current ripple as a function of current


def RippleFunction(I0):
    RipMax = 3.0
    RipFS = 1.0
    FS = 12.5e3
    dI = (RipMax - RipFS) / FS**2 * (-(I0 - FS) * (I0 + FS)) + RipFS
    return dI


if False:
    dRdB = []  # change in Target position with respect to B
    fB = 0.01
    for i in range(len(Bn)):
        B = BfieldTF(B0=Bn[i])
        Bv = BfieldVF(B0=0.00000)
        T = Trajectory(Vessel, B, Bv, v0=Vinjection, E0=Energy)
        T.LineColor = CMAP(1.0 * i / len(Bn))
        T.LineWidth = 2.0
        TrajectoryList.append(T)

        B1 = BfieldTF(B0=Bn[i] * (1 + fB))  # Bfield + fractional change
        T1 = Trajectory(Vessel, B1, Bv, v0=Vinjection, T0=Energy)
        dRdB.append(T1.Target.Distance(T) / (fB * 100.0))
    print(dRdB)

# ------------------------------------------------------------------------------
# Calculate Target Error trajectories given % Magnet Ripple
if True:
    #	fB = 0.02
    # Centroid Trajectory
    I0 = []
    for i in range(len(Bn)):
        B = BfieldTF(B0=Bn[i])
        Bv = BfieldVF(B0=0.00000)
        T = Trajectory(Vessel, B, Bv, v0=Vinjection, T0=Energy)
        T.LineColor = CMAP(1.0 * i / len(Bn))
        T.LineWidth = 1.0
        T.LineStyle = ':'
        I0.append(CalculateI0(Bn[i]))
        TrajectoryList.append(T)

# Error Trajectories
    RError = [np.zeros(3), np.zeros(3)]
    DeltaR = []
    for i in range(len(Bn)):
        I1 = CalculateI0(Bn[i])
        dI = RippleFunction(I1) / 2.0
#		fB = CalculateB0(dI)
        fB = CalculateB0(dI + 0.0005 * 12.5e3)
        fB1 = [-fB, fB]
        print(dI, fB)
        for j in range(len(fB1)):
            fB1 = [-fB, fB]
#			B1 = BfieldTF(B0=Bn[i]*(1+fB1[j]) ) # B * (1 + epsilon)
            B1 = BfieldTF(B0=Bn[i] + fB1[j])  # B + epsilon
            T1 = Trajectory(Vessel, B1, Bv, v0=Vinjection,
                            T0=Energy, dS=1e-4, Nmax=100000)
            RError[j] = T1.target.XYZ
            T1.LineStyle = '-'
            T1.LineWidth = 1.0
            T1.LineColor = CMAP(1.0 * i / len(Bn))
#			TrajectoryList.append(T1)
        DeltaR.append(norm(RError[1] - RError[0]))

    # Save Target parameters
#	T.Target.SaveTargetParameters(TFCurrent=In[i],Path=OutputPath+'geometry/')

    # append lists of Target Quantities
#	AngleComponents.append([T.Target.VAngle,T.Target.HAngle])
#	Coordinates.append([T.Target.R,T.Target.Z,T.Target.Phi])
#	Parameters.append(T.Target.GetDetectionParameters())

# ===============================================================================
# Plot Resuts
# ===============================================================================
if True:
    # Plot 3D results

    for i in range(len(TrajectoryList)):
        TrajectoryList[i].Plot3D(ax)
#		TrajectoryList[i].Target.Plot3D(ax);
# TrajectoryList[-1].Limits3D(ax);

# ------------------------------------------------------------------------------
# Construct Legend
    Leg = []
    for i in range(len(Bn)):
        #		Leg.append('B = %0.3fT' % TrajectoryList[i].BFieldTF.B0)
        Leg.append('B = %0.2f kA' % ((TrajectoryList[i].BFieldTF.I0) / 1000.0))

# ------------------------------------------------------------------------------
# Plot 2D projections of trajectories (Poloidal View)
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
    plt.title(r'Poloidal Projection ($\alpha$ = %0.1f$^o$, $\beta$ = %0.1f$^o$)' % (
        alpha0, beta0))
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
    plt.title(r'Midplane Projection ($\alpha$ = %0.1f$^o$, $\beta$ = %0.1f$^o$)' % (
        alpha0, beta0))
    ax = plt.subplot(1, 2, 2)
    ax.legend(Leg, bbox_to_anchor=(1.28, 1.0))

#	plt.legend(('B = 0.05','B = 0.10','B = 0.15','B = 0.20','B = 0.25','B = 0.30')

# ------------------------------------------------------------------------------
# Plot Error vs I0
plt.figure()
plt.plot(np.array(I0) * 1e-3, np.array(DeltaR) * 1e3)
plt.xlabel('Current [kA]')
plt.ylabel('peak-peak Error [mm]')
plt.title('Beam centroid error due to TF current error')

# ===============================================================================
# Save Angular and Detection Quantities
# ===============================================================================

if False:
    np.savetxt(OutputPath + 'geometry/TargetAngle_Vert_Horiz.dat', AngleComponents)
    np.savetxt(OutputPath + 'geometry/TargetCoordinates.dat', Coordinates)
    Header0 = '(0) I0 [A], (1) B0 [T], (2) X [m] , (3) Y [m], (4) Z [m], (5) incident angle [rad], (6) Detection Angle [rad], (7) optical path length [m] , (8) Detection Angle [rad], (9) Detection Angle [deg], (10) Detector Eff'
    np.savetxt(OutputPath + 'geometry/DetectionParameters.dat',
            (np.array(Parameters)), header=Header0)

if False:
    FigName = 'TrajectoryProjections_alpha%2.2f_beta%2.2f_' % (
        alpha0, beta0)  # + B.Method
    FigPath = '../output/plots/'
    TrajectoryList[-1].Target.SaveTargetParameters(
        Path=FigPath + 'Test_alpha%2.2f_beta%2.2f_UpDown' % (alpha0, beta0))
    plt.savefig(FigPath + FigName + '_UpDown.pdf')
    plt.savefig(FigPath + FigName + '_UpDown.png')
    print('File saved: ' + FigName)


plt.show()
