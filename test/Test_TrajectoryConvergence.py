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

Rb = [0.1, 0.1, 100.0, 100.0]
Zb = [-100.0, 100.0, 100.0, -100.0]

Vessel = Boundary(Rb, Zb)
# Vessel.Plot2D(0)

if False:
    B = Bfieldc(B0=0.1)
    Bv = Bfieldc(B0=0.0001)
    d0 = 10.0
    dS = logspace(-5, -2, 15)
    dr = []
    T = []
    for i in range(len(dS)):
        #		T.append(Trajectory(Vessel,B,r0=[20.0,0.0,0.0],v0=[0.0,0.0,1.0],dS=dS[i],Nmax=round(d0/dS[i])) )
        T.append(Trajectory(Vessel, B, Bv, r0=[20.0, 0.0, 0.0], v0=[
                 0.0, 0.0, 1.0], dS=dS[i], Nmax=round(d0 / dS[i])))
        RL = (T[-1].m0 * T[-1].v0) / (T[-1].q0 * B.B0)
        R = T[-1].r[-1] - np.array([20.0 - RL, 0, 0.0])
        # -RL)/d0 ) #/T.s[-1]*d0 - RL)
        dr.append(np.sqrt(R[0]**2 + R[1]**2 + R[2]**2) * (d0 / T[-1].s[-1]) - RL)
        # T.Plot2D()
    plt.figure(1)
    plt.loglog(dS, dr, '.')
    plt.xlabel(r'Step Size $\Delta$S [m]')
    plt.ylabel(r'$\Delta R/S$')
    plt.title(r'Error / Arc length')

    plt.figure(2)
    plt.loglog(dS, dr / RL, '.')
    plt.xlabel(r'Step Size $\Delta$S [m]')
    plt.ylabel(r'${\Delta R}/{R S}$ [1/m]')
    plt.title(r'Normalized Bending Radius Error per Arc Length')

# Test Varying Total Distance S
if True:
    BTF = [0.15, 0.1, 0.05]
    for i in range(len(BTF)):
        B = Bfieldc(B0=BTF[i])
        Bv = Bfieldc(B0=0.0)
        T = Trajectory(Vessel, B, Bv, r0=[20.0, 0.0, 0.0], v0=[
                       0.0, 1.0, 1.0], Nmax=10000)
        # T.Plot2D()

        x = []
        y = []
        z = []
        S = []
        r = []
        R = []
        rN = []
        # RL = (T.m0*T.v0) / (T.q0*B.B0)  # mV/qB
        RL = T.m0 / T.c0 * T.beta[-1] / (B.B0)  # mV/qB
        R0 = np.array([20.0 - RL, 0, 0.0])

        for i in range(len(T.r)):
            R.append(T.r[i] - R0)

        for i in range(len(T.r)):
            x.append(R[i][0])
            y.append(R[i][1])
            z.append(R[i][2])
            S.append(T.s[i])
            r.append(np.sqrt(x[-1]**2 + y[-1]**2 + z[-1]**2))
            rN.append((r[-1] - RL) / T.s[i])

        plt.figure(1)
        plt.plot(x, z)
        plt.figure(2)
        plt.plot(S, abs(1 - r / RL))
        plt.xlabel('S Coordinate [m]')
        plt.ylabel(r'Error $\epsilon $')
        plt.title(r'$\epsilon = \Delta $r/Rc')
#		plt.figure(3); plt.loglog(S,abs(1-r/RL));
#		plt.xlabel('S Coordinate [m]'); plt.ylabel(r'$\epsilon $'); plt.title(r'$\epsilon = \Delta $r/Rc')
        plt.figure(4)
        plt.plot(S, np.array(rN))
        plt.xlabel('S Coordinate [m]')
        plt.ylabel(r'Normalized Error $\epsilon_N $')
        plt.title(r'$\epsilon_N = \Delta$ r/S')
        plt.figure(5)
        plt.semilogx(S, np.array(rN))
        plt.xlabel('S Coordinate [m]')
        plt.ylabel(r'Normalized Error $\epsilon_N $')
        plt.title(r'$\epsilon_N = \Delta$ r/S')

        plt.figure(6)
        plt.subplot(2, 1, 1)
        plt.plot(S, abs(1 - r / RL))
        plt.xlabel('S Coordinate [m]')
        plt.ylabel(r'Error $\epsilon_r$')
        plt.title(r'Radius of Curvature Error $\epsilon_r = \Delta$r/Rc')
        plt.xlim(0, 100)
        plt.subplot(2, 1, 2)
        plt.semilogx(S, np.array(rN))
        plt.xlabel('S Coordinate [m]')
        plt.ylabel(r'Error $\epsilon_N $')
        plt.title(
            r'Transverse Error Normalized to Trajectory Length $\epsilon_N = \Delta$ r/S')
        plt.xlim(0.1, 100)
    plt.subplot(2, 1, 1)
    plt.legend((r'B$_{TF}$ = 0.15 T', r'B$_{TF}$ = 0.10 T',
               r'B$_{TF}$ = 0.05 T'), loc=2)


plt.show()
