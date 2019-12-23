
import sys
sys.path.append('../lib/')
from BeamDynamicsTools import *
import numpy as np
import pylab as pl

L0 = 0.1
L1 = 10.0

Rb = [L0, L0, L1, L1, L0]
Zb = [0.0, L1, L1, -L1, -L1]

Vessel = Boundary(Rb, Zb)
Vessel.Plot2D(0)

# R = m v / q B -> v = np.sqrt(2 E / m) -> np.sqrt( 2 m E) / q B

#DATA = loadtxt('CmodCoordinatesRZ.txt')
#Rb=[]; Zb=[];
# for i in range(len(DATA[:,0])):
#	Rb.append(DATA[i,0])
#	Zb.append(DATA[i,1])

Vessel = Boundary(Rb, Zb)
Vessel.Plot2D(0)

# R = np.sqrt(2 M E)*(c/B)
R0 = 1.0
Angle = np.pi / 2.0
BR = np.sqrt(2.0 * (2.0 * 1.67262158e-27) *
          (0.9e6 * 1.602e-19)) / (1.60217646e-19 * R0)
Bv = BfieldVF(B0=0.00000)

#DeltaS = R0*Angle
if True:
    #	B = BfieldTF(B0=0.2)
    B = Bfieldc(B0=BR)
    T = Trajectory(Vessel, B, Bv, r0=[
                   10.0, 0.0, 1.0], v0=[-1.0, 0.0, 0.0], Nmax=100)
    plot(Rb, Zb)
    T.Plot2D()
#	T.PlotB()
    pl.xlabel(r'X')
    pl.ylabel(r'Y')
    pl.title(r'?')
    xlim(-L0, L1 * 1.1)
    ylim(-L1 * 1.1, L1 * 1.1)
    Y = []
    for i in range(len(T.B)):
        #		Y.append(norm( cross(T.v[i],T.B[i]) ) )
        Yi = np.cross(T.v[i], T.B[i])
        Y.append(Yi[2])

#	pl.figure(); pl.plot(Y)

if False:
    S0 = np.matrix([
        [0.577100, 0.398000, 0.000000, 0.000000, 0.000000, 0.000000],
        [0.398000, 171.8262, 0.000000, 0.000000, 0.000000, 0.000000],
        [0.000000, 0.000000, 0.343900, -0.27150, 0.000000, 0.000000],
        [0.000000, 0.000000, -0.27150, 238.3722, 0.000000, 0.000000],
        [0.000000, 0.000000, 0.000000, 0.000000, 1.297156, 2.343722],
        [0.000000, 0.000000, 0.000000, 0.000000, 2.343722, 134.9344]], float)

    S0 = np.matrix(loadtxt('data/SigmaFinal/SigmaInj.txt'))

    Beam = beam(T, S0, Target=False)

    Beam.Trace(Target=False)

    pl.figure
    STR = r'1m Radius 90$^o$ Bend in Uniform B$_x$- Field'
    Ei = Ellipse(Beam.Sigma[0])
    Ei.PlotALL()
    Em1 = Ellipse(Beam.Sigma[int(len(Beam.Sigma) * 0.33)])
    Em1.PlotALL()
    Em2 = Ellipse(Beam.Sigma[int(len(Beam.Sigma) * 0.66)])
    Em2.PlotALL()
    Ef = Ellipse(Beam.Sigma[-1])
    Ef.PlotALL(Title=STR)

#	savetxt('data/SigmaFinal/SigmaBend90.txt',Beam.Sigma[-1])

pl.show()
