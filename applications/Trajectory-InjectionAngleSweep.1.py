import numpy as np
import matplotlib.pyplot as plt
import sys
import os
sys.path.append(os.path.join('../'))

from lib.BeamDynamicsTools.Boundary import Boundary
from lib.BeamDynamicsTools.Bfield import Bfield, BfieldTF, BfieldVF
from lib.BeamDynamicsTools.Trajectory import Trajectory
from lib.BeamDynamicsTools.Beam import Beam


if __name__ == "__main__":
    # Define np.array of injection angles
    #  (x,y,z) = (1.798m, -0.052m, 0.243m)
    #  alpha = 12.6 degrees (X-Z plane)
    #  beta = 8.0 degrees (X-Y plane)
    #  VInjection
    alpha = np.linspace(8, 12, 4) / 180.0 * np.pi
    beta = 8.0 / 180.0 * np.pi
    Rinjection = [1.798, 0.052, 0.243]
    Vinjection = []
    for i, alph in enumerate(alpha):
        v0 = -np.cos(alph) * np.cos(beta)
        v1 = np.cos(alph) * np.sin(beta)
        v2 = -np.sin(alph)
        Vinjection.append([v0, v1, v2])

    # ------------------------------------------------------------------------------
    # Import poloidal Boundary points
    Rb = np.loadtxt('../data/CmodCoordinatesRZ.dat', usecols=[0])
    Zb = np.loadtxt('../data/CmodCoordinatesRZ.dat', usecols=[1])

    # ------------------------------------------------------------------------------
    # Generate vessel Boundary
    Vessel = Boundary(Rb, Zb)

    # ------------------------------------------------------------------------------
    # 3D plot of vessel Boundary
    ax = Vessel.Figure3D(1)
    Vessel.Plot3D(ax)

    # ------------------------------------------------------------------------------
    # Inputs for B-field settings
    In = np.array([0.0, 1600.0, 3120, 4450.0])
    Bn = np.array([0.0, 0.05818182, 0.11345455, 0.16181818])

    # ===============================================================================
    # Perform Trajectory calculation for Trajectory Sweep
    # ===============================================================================

    AngleComponents = []
    Coordinates = []
    Parameters = []
    TrajectoryList = []
    OutputPath = './tmp/'
    Color = ['b', 'g', 'r', 'c']
    for j, alph in enumerate(alpha):
        for i, B0 in enumerate(Bn):
            B = BfieldTF(B0=B0)
            Bv = BfieldVF(B0=0.00000)
            T = Trajectory(Vessel, B, Bv, v0=Vinjection[j], Method='LeapFrog')
            T.LineColor = Color[j]
            T.Plot3D(ax)
            T.PlotB(FIG=2)
            T.PlotV(FIG=3)
            TrajectoryList.append(T)

            np.savetxt(
                "./tmp/Traject_{:d}_{:d}.dat".format(j, i), np.array(T.r))
            AngleComponents.append([T.target.VAngle, T.target.HAngle])
            Coordinates.append([T.target.R, T.target.Z, T.target.Phi])
            Parameters.append(T.target.GetDetectionParameters())

    T.Limits3D(ax)

    plt.figure(figsize=(20, 8))
    for T in TrajectoryList:
        plt.subplot(1, 2, 1)
        T.Plot2D('poloidal')
    for T in TrajectoryList:
        plt.subplot(1, 2, 2)
        T.Plot2D('top')

    # ------------------------------------------------------------------------------
    # Save Angular and Detection Quantities
    np.savetxt(OutputPath + 'TargetAngle_Vert_Horiz.dat', AngleComponents)
    np.savetxt(OutputPath + 'TargetCoordinates.dat', Coordinates)
    Header0 = '(0) I0 [A], (1) B0 [T], (2) X [m] , (3) Y [m], (4) Z [m], (5) incident angle [rad], (6) Detection Angle [rad], (7) optical path length [m] , (8) Detection Angle [rad], (9) Detection Angle [deg], (10) Detector Eff'
    np.savetxt(OutputPath + 'DetectionParameters.dat',
               (np.array(Parameters)), header=Header0)

    plt.show()
