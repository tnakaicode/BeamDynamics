# Trajectory.py
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as cnt
import os
import sys
import math
import timeit
from mpl_toolkits.mplot3d import Axes3D
from numpy.linalg import norm

sys.path.append(os.path.join('../'))
from .Ellipse import *
from .AngleCorrection import *
from .Target import *


# ======= Default injection geometry ==================
# (x,y,z) = (1.798m, -0.052m, 0.243m)
#  alpha = 12.6 degrees (X-Z plane)
#  beta = 8.0 degrees (X-Y plane)
alpha = 12.6 / 180.0 * np.pi
beta = 8.0 / 180.0 * np.pi

Rinjection = [1.798, -0.052, 0.243]
Vinjection = [-np.cos(alpha) * np.cos(beta), np.cos(alpha)
              * np.sin(beta), -np.sin(alpha)]

Mass0 = 2.0 * (938.272e6)
dLB = 2.0e-3
# scale length for B gradient
#Vinjection = [-1,0,0]
# ====== \Default injection geometry ==================


class Trajectory(object):

    def __init__(self, Vessel, B, Bv, dS=1e-3, r0=Rinjection, v0=Vinjection, a0=[0.0, 0.0, 0.0], M0=Mass0, T0=0.9e6, I0=1e-3, Freq=425e6, Nmax=5000, Smin=1.1, Smax=5.0, Method='Relativistic'):
        self.start = timeit.default_timer()

        self.fig_3d = plt.figure(1)
        self.axs_3d = Axes3D(self.fig_3d)

        self.fig_2d, self.axs_2d = plt.subplots(1, 2, figsize=(15, 8))
        self.axs_pol = self.axs_2d[0]
        self.axs_pol.set_aspect('auto')
        self.axs_pol.xaxis.grid()
        self.axs_pol.yaxis.grid()
        self.axs_top = self.axs_2d[1]
        self.axs_top.set_aspect('auto')
        self.axs_top.xaxis.grid()
        self.axs_top.yaxis.grid()

        self.FigureBV()

        # Method == 'Relativistic'
        # Method == 'LeapFrog'
        # Method == 'Euler'

        # B = Magnetic Field [T] (BfieldTF class)
        # Vessel = Defines wall (Boundary class)
        # M0 = Rest Mass [MeV/c^2]
        # T0 = kinetic energy beam [eV]
        # r  = position vector [x, y, z]
        # v  = velocity vector [Vx, Vy, Vz]
        # a  = acceleration vector [ax, ay, az]
        # Proton Mass: 938.272e6
        # Electron Charge: 1.60217646e-19

        # Particle and beam constants
        self.c0 = cnt.c
        self.A0 = M0 / 938.272e6
        self.q0 = 1.60217646e-19
        self.qm = self.q0 / (M0 * self.q0 / self.c0**2)
        self.m0 = M0
        self.I0 = I0
        self.Frequency = Freq
        self.T0 = T0

        # Vessel Boundary
        self.Vessel = Vessel

        # Magnetic coil sets
        self.BFieldTF = B
        self.BFieldVF = Bv

        # Beam
        self.r = [np.array(r0)]
        self.gamma = 1.0 + T0 / M0
        self.beta = np.sqrt(1.0 - 1.0 / self.gamma**2)
        self.Beta = [self.beta * np.array(v0) / norm(v0)]
#		self.v0 = np.sqrt(2.0*T0*1.602e-16/(A0*1.67e-27))
#		self.v0 = np.sqrt(2.0*T0*self.q0/(self.m0))
#		self.v0 = np.sqrt(2.0*T0*self.q0/(self.m0))
        self.v0 = self.beta * self.c0
        self.v = [self.v0 * np.array(v0) / norm(v0)]
        self.beta = [norm(self.Beta[-1])]
        self.gamma = [1.0 / (1.0 - self.beta[-1]**2)]
        self.a = [np.array(a0)]
        self.F = [0.0]
        self.B = [np.array(B.local(r0))]
        self.s = [0.0]
        self.dt = dS / self.v0
        self.dS = [dS]

        # Gradient and curvature attributes
        self.k = [0.0]
        self.Rc = []
        self.gradB = [0.0]
        self.gradBk = [0.0]
        self.gradBn = [0.0]
        self.gradBx = [0.0]
        self.gradBy = [0.0]

        # Plotting attributes
        self.LineColor = 'r'
        self.LineWidth = 2
        self.LineStyle = '-'

        # num Limit
        self.Nmax = Nmax
        self.Smax = Smax
        self.Smin = Smin

        if Method == 'Relativistic':
            self.Method_Relativistic()
        elif Method == 'LeapFrog':
            self.Method_LeapFrog()
        elif Method == 'Euler':
            self.Method_Euler()
        elif Method == None:
            pass
        else:
            self.Method_Relativistic()

    def init_condition(self, Vessel, B, Bv, dS=1e-3, r0=Rinjection, v0=Vinjection, a0=[0.0, 0.0, 0.0], M0=Mass0, T0=0.9e6, I0=1e-3, Freq=425e6, Nmax=5000, Smin=1.1, Smax=5.0, Method='Relativistic'):
        self.start = timeit.default_timer()

        # B = Magnetic Field [T] (BfieldTF class)
        # Vessel = Defines wall (Boundary class)
        # M0 = Rest Mass [MeV/c^2]
        # T0 = kinetic energy beam [eV]
        # r  = position vector [x, y, z]
        # v  = velocity vector [Vx, Vy, Vz]
        # a  = acceleration vector [ax, ay, az]
        # 938.272e6 Proton Mass
        # Electron Charge

        # Particle and beam constants
        self.c0 = cnt.c
        self.A0 = M0 / 938.272e6
        self.q0 = 1.60217646e-19
        self.qm = self.q0 / (M0 * self.q0 / self.c0**2)
        self.m0 = M0
        self.I0 = I0
        self.Frequency = Freq
        self.T0 = T0

        # Vessel Boundary
        self.Vessel = Vessel

        # Magnetic coil sets
        self.BFieldTF = B
        self.BFieldVF = Bv

        # Beam
        self.r = [np.array(r0)]
        self.gamma = 1.0 + T0 / M0
        self.beta = np.sqrt(1.0 - 1.0 / self.gamma**2)
        self.Beta = [self.beta * np.array(v0) / norm(v0)]
#		self.v0 = np.sqrt(2.0*T0*1.602e-16/(A0*1.67e-27))
#		self.v0 = np.sqrt(2.0*T0*self.q0/(self.m0))
#		self.v0 = np.sqrt(2.0*T0*self.q0/(self.m0))
        self.v0 = self.beta * self.c0
        self.v = [self.v0 * np.array(v0) / norm(v0)]
        self.beta = [norm(self.Beta[-1])]
        self.gamma = [1.0 / (1.0 - self.beta[-1]**2)]
        self.a = [np.array(a0)]
        self.F = [0.0]
        self.B = [np.array(B.local(r0))]
        self.s = [0.0]
        self.dt = dS / self.v0
        self.dS = [dS]

        # Gradient and curvature attributes
        self.k = [0.0]
        self.Rc = []
        self.gradB = [0.0]
        self.gradBk = [0.0]
        self.gradBn = [0.0]
        self.gradBx = [0.0]
        self.gradBy = [0.0]

        # Plotting attributes
        self.LineColor = 'r'
        self.LineWidth = 2
        self.LineStyle = '-'

        # num Limit
        self.Nmax = Nmax
        self.Smax = Smax
        self.Smin = Smin

        if Method == 'Relativistic':
            self.Method_Relativistic()
        elif Method == 'LeapFrog':
            self.Method_LeapFrog()
        elif Method == 'Euler':
            self.Method_Euler()
        elif Method == None:
            pass
        else:
            self.Method_Relativistic()


# ===============================================================================
# Euler Integration:
# ===============================================================================


    def Method_Euler(self):
        c1 = True
        c2 = True
        i = 0
        IN = False
        NormalV = np.zeros(3)
        TangentV = np.zeros(3)
        IncidentV = np.zeros(3)
        RT = np.zeros(3)

        while (c1 or c2) and (i < self.Nmax and self.s[-1] < self.Smax):
            self.r.append(self.r[-1] + self.v[-1] * self.dt)
            self.s.append(self.s[-1] + self.dS[-1])
            self.B.append(self.BFieldTF.local(
                self.r[-1]) + self.BFieldVF.local(self.r[-1]))
            self.a.append(self.qm * np.cross(self.v[-1], self.B[-1]))
            self.v.append(self.v[-1] + (self.a[-1]) * self.dt)
            self.dS.append(self.s[-1] - self.s[-2])
# ------------------------------------------------------------------------------
            # Normalized Relativistic Parameters
            self.Beta.append(self.v[-1] / self.c0)
            self.beta.append(norm(self.Beta[-1]))
            self.gamma.append(1.0 / (1.0 - self.beta[-1]**2))
# ------------------------------------------------------------------------------
            # Check to see if beam np.crosses boundary
            IN = True
            c3 = self.s[-1] > self.Smin
            c4 = self.Vessel.InBoundary(self.r[-1])
            c5 = self.s[-1] < self.Smax
            if c3:
                if c4:
                    IN, NormalV, TangentV, IncidentV, RT, Xpol = self.Vessel.Xboundary(
                        self.r[-2], self.r[-1])
# ------------------------------------------------------------------------------
            # record curvature and bending radius
            #self.k.append(qm * cross(self.v[-1],self.B[-1])/self.v0**2)
            self.k.append(norm(self.a[-1] / self.v0**2))
            self.Rc.append(1.0 / self.k[-1])
# ------------------------------------------------------------------------------
            # B Record Gradients
            vecR = -1.0 * (self.a[-1]) / norm(self.a[-1])
            vecB = self.B[-1] / norm(self.B[-1])
            Br2 = norm(self.BFieldTF.local(self.r[-1] + vecR * dLB))
            Br1 = norm(self.BFieldTF.local(self.r[-1] - vecR * dLB))
            Bb2 = norm(self.BFieldTF.local(self.r[-1] + vecB * dLB))
            Bb1 = norm(self.BFieldTF.local(self.r[-1] - vecB * dLB))
            # ( np.array( [(Br2-Br1)/(2.0*dLB) , (Bb2-Bb1)/(2.0*dLB)] ) )
            self.gradB.append((Br2 - Br1) / (2.0 * dLB))
            # (qm/(self.gamma[-1]*self.beta[-1]*c0**2))
            self.gradBk.append(self.gradB[-1] * self.qm / (self.c0 * self.v0))
            self.gradBn.append(-1.0 * self.Rc[-1] /
                               norm(self.B[-1]) * (Br2 - Br1) / (2.0 * dLB))
            self.gradBx.append((Br2 - Br1) / (2.0 * dLB))
            self.gradBy.append((Bb2 - Bb1) / (2.0 * dLB))
# ------------------------------------------------------------------------------
            # Conditional statements for continuing iteration
            c1 = IN
            c2 = self.s[-1] < self.Smin
            i = i + 1
# ------------------------------------------------------------------------------
        # Calculate Basis Matrices
        self.BeamBasis()
        self.stop = timeit.default_timer()
        self.RunTime = self.stop - self.start
        print('trajectory complete, S = %0.3f m, B0 = %0.4f T, B0 = %0.4f T, RunTime = %0.1f s' % (
            self.s[-1], self.BFieldTF.B0, self.BFieldVF.B0, self.RunTime))
# ------------------------------------------------------------------------------
        # Define Target
        if i < self.Nmax - 1 and self.s[-1] <= self.Smax:
            self.target = Target(NormalV, TangentV, IncidentV,
                                 self.BFieldTFBFieldTF, self.BFieldTFBFieldVF, RT, Xpol)
            self.target.SigmaBasis = self.BasisM6[-1]
# ----------------------------------------------------------------------------
        # If no boundary was reached assume normal incidence
        else:
            NormalV = np.array(self.BasisM3[-1][:, 2]).flatten()
            TangentV = np.array(self.BasisM3[-1][:, 1]).flatten()
            IncidentV = np.array(self.BasisM3[-1][:, 2]).flatten()
            RT = self.r[-1]
            self.target = Target(NormalV, TangentV, IncidentV,
                                 self.BFieldTF, self.BFieldVF, RT, Xpol)
        print('Beam Coordinates Complete')
# ===============================================================================
# End
# Euler Integration:
# ===============================================================================

# ===============================================================================
# Leapfrog Integration:
# ===============================================================================
    def Method_LeapFrog(self):
        c1 = True
        c2 = True
        i = 0
        IN = False
        NormalV = np.zeros(3)
        TangentV = np.zeros(3)
        IncidentV = np.zeros(3)
        RT = np.zeros(3)

        while (c1 or c2) and (i < self.Nmax and self.s[-1] < self.Smax):
            self.r.append(self.r[-1] + self.v[-1] *
                          self.dt + 0.5 * self.a[-1] * self.dt**2)
            self.s.append(self.s[-1] + self.dS[-1])
            self.B.append(self.BFieldTF.local(
                self.r[-1]) + self.BFieldVF.local(self.r[-1]))
            self.a.append(self.qm * np.cross(self.v[-1], self.B[-1]))
            self.v.append(self.v[-1] + 0.5 *
                          (self.a[-1] + self.a[-2]) * self.dt)
            self.dS.append(self.s[-1] - self.s[-2])
# ------------------------------------------------------------------------------
            # Normalized Relativistic Parameters
            self.Beta.append(self.v[-1] / self.c0)
            self.beta.append(norm(self.Beta[-1]))
            self.gamma.append(1.0 / (1.0 - self.beta[-1]**2))
# ------------------------------------------------------------------------------
            # Check to see if beam crosses boundary
            IN = True
            c3 = self.s[-1] > self.Smin
            c4 = self.Vessel.InBoundary(self.r[-1])
            c5 = self.s[-1] < self.Smax
            if c3:
                if (not c4):
                    IN, NormalV, TangentV, IncidentV, RT, Xpol = self.Vessel.Xboundary(
                        self.r[-2], self.r[-1])
# ------------------------------------------------------------------------------
            # record bending radius
            # self.k.append(qm * cross(self.v[-1],self.B[-1])/self.v0**2)
            self.k.append(norm(self.a[-1] / self.v0**2))
            self.Rc.append(1.0 / self.k[-1])
# ------------------------------------------------------------------------------
            # B Record Gradients
            vecR = -1.0 * (self.a[-1]) / norm(self.a[-1])
            vecB = self.B[-1] / norm(self.B[-1])
            Br2 = norm(self.BFieldTF.local(self.r[-1] + vecR * dLB))
            Br1 = norm(self.BFieldTF.local(self.r[-1] - vecR * dLB))
            Bb2 = norm(self.BFieldTF.local(self.r[-1] + vecB * dLB))
            Bb1 = norm(self.BFieldTF.local(self.r[-1] - vecB * dLB))
            # ( np.array( [(Br2-Br1)/(2.0*dLB) , (Bb2-Bb1)/(2.0*dLB)] ) )
            self.gradB.append((Br2 - Br1) / (2.0 * dLB))
            # (qm/(self.gamma[-1]*self.beta[-1]*c0**2))
            self.gradBk.append(self.gradB[-1] * self.qm / (self.c0 * self.v0))
            self.gradBn.append(-1.0 * self.Rc[-1] /
                               norm(self.B[-1]) * (Br2 - Br1) / (2.0 * dLB))
            self.gradBx.append((Br2 - Br1) / (2.0 * dLB))
            self.gradBy.append((Bb2 - Bb1) / (2.0 * dLB))
# ------------------------------------------------------------------------------
            # Conditional statements for continuing iteration
            c1 = IN
            c2 = self.s[-1] < self.Smin
            i = i + 1
# ------------------------------------------------------------------------------
        # Calculate Basis Matrices
        self.BeamBasis()
        self.stop = timeit.default_timer()
        self.RunTime = self.stop - self.start
        print('trajectory complete, S = %0.3f m, B0 = %0.4f T, B0 = %0.4f T, RunTime = %0.1f s' % (
            self.s[-1], self.BFieldTF.B0, self.BFieldVF.B0, self.RunTime))
# ------------------------------------------------------------------------------
        # Define Target
        if i < self.Nmax - 1 and self.s[-1] <= self.Smax:
            self.target = Target(
                NormalV, TangentV, IncidentV, self.BFieldTF, self.BFieldVF, RT, Xpol)
            self.target.SigmaBasis = self.BasisM6[-1]
# ------------------------------------------------------------------------------
        # If no boundary was reached assume normal incidence
        else:
            NormalV = np.array(self.BasisM3[-1][:, 2]).flatten()
            TangentV = np.array(self.BasisM3[-1][:, 1]).flatten()
            IncidentV = np.array(self.BasisM3[-1][:, 2]).flatten()
            RT = self.r[-1]
            self.target = Target(
                NormalV, TangentV, IncidentV, self.BFieldTF, self.BFieldVF, RT, Xpol)
        print('Beam Coordinates Complete')
# ===============================================================================
# End
# Leapfrog Integration:
# ===============================================================================

# ===============================================================================
# Relativistic Euler Integration:
# ===============================================================================
    def Method_Relativistic(self):
        c1 = True
        c2 = True
        i = 0
        IN = False
        NormalV = np.zeros(3)
        TangentV = np.zeros(3)
        IncidentV = np.zeros(3)
        RT = np.zeros(3)

        while (c1 or c2) and (i < self.Nmax and self.s[-1] < self.Smax):
            self.r.append(self.r[-1] + self.v[-1] * self.dt)
            self.s.append(self.s[-1] + self.dS[-1])
            self.B.append(self.BFieldTF.local(
                self.r[-1]) + self.BFieldVF.local(self.r[-1]))
            self.F.append(self.q0 * np.cross(self.v[-1], self.B[-1]))
            self.a.append(self.c0**2 / (self.gamma[-1] * self.m0 * self.q0) * (
                self.F[-1] - (np.dot(self.v[-1], self.F[-1]) * self.v[-1] / self.c0**2)))
            self.v.append(self.v[-1] + self.a[-1] * self.dt)
            self.dS.append(self.s[-1] - self.s[-2])
# ------------------------------------------------------------------------------
            # Normalized Relativistic Parameters
            self.Beta.append(self.v[-1] / self.c0)
            self.beta.append(norm(self.Beta[-1]))
            self.gamma.append(1.0 / (1.0 - self.beta[-1]**2))

            # Check to see if beam crosses boundary
            IN = True
            c3 = self.s[-1] > self.Smin
            c4 = self.Vessel.InBoundary(self.r[-1])
            c5 = self.s[-1] < self.Smax
            if c3:
                if (not c4):
                    IN, NormalV, TangentV, IncidentV, RT, Xpol = self.Vessel.Xboundary(
                        self.r[-2], self.r[-1])
# ------------------------------------------------------------------------------
            # record curvature and bending radius
            self.k.append(norm(self.a[-1] / self.v0**2))
            self.Rc.append(1.0 / self.k[-1])

# ------------------------------------------------------------------------------
            # B Record Gradients
            vecR = -1.0 * (self.a[-1]) / norm(self.a[-1])
            vecB = self.B[-1] / norm(self.B[-1])
            Br2 = norm(self.BFieldTF.local(self.r[-1] + vecR * dLB))
            Br1 = norm(self.BFieldTF.local(self.r[-1] - vecR * dLB))
            Bb2 = norm(self.BFieldTF.local(self.r[-1] + vecB * dLB))
            Bb1 = norm(self.BFieldTF.local(self.r[-1] - vecB * dLB))
            # ( np.array( [(Br2-Br1)/(2.0*dLB) , (Bb2-Bb1)/(2.0*dLB)] ) )
            self.gradB.append((Br2 - Br1) / (2.0 * dLB))
            # (qm/(self.gamma[-1]*self.beta[-1]*c0**2))
            self.gradBk.append(self.gradB[-1] * self.qm / (self.c0 * self.v0))
            self.gradBn.append(-1.0 * self.Rc[-1] /
                               norm(self.B[-1]) * (Br2 - Br1) / (2.0 * dLB))
            self.gradBx.append((Br2 - Br1) / (2.0 * dLB))
            self.gradBy.append((Bb2 - Bb1) / (2.0 * dLB))

# ------------------------------------------------------------------------------
            # Conditional statements for continuing iteration
            c1 = IN
            c2 = self.s[-1] < self.Smin
            i = i + 1

# ------------------------------------------------------------------------------
        # Calculate Basis Matrices
        self.BeamBasis()
        self.stop = timeit.default_timer()
        self.RunTime = self.stop - self.start
        print('trajectory complete, S = %0.3f m, B0 = %0.4f T, B0 = %0.4f T, RunTime = %0.1f s' % (
            self.s[-1], self.BFieldTF.B0, self.BFieldVF.B0, self.RunTime))
# ------------------------------------------------------------------------------
        # Define Target
        if i < self.Nmax - 1 and self.s[-1] <= self.Smax:
            self.target = Target(NormalV, TangentV, IncidentV,
                                 self.BFieldTF, self.BFieldVF, RT, Xpol)
            self.target.SigmaBasis = self.BasisM6[-1]
            #self.target.TargetBasis = self.BasisM6[-1]
# ------------------------------------------------------------------------------
        # If no boundary was reached assume normal incidence
        else:
            NormalV = np.array(self.BasisM3[-1][:, 2]).flatten()
            TangentV = np.array(self.BasisM3[-1][:, 1]).flatten()
            IncidentV = np.array(self.BasisM3[-1][:, 2]).flatten()
            RT = self.r[-1]
            self.target = Target(NormalV, TangentV, IncidentV,
                                 BFieldTF, BFieldVF, RT, Xpol)
        print('Beam Coordinates Complete')
# ===============================================================================
# End
# Relativistic Euler Integration:
# ===============================================================================

    # Calculate 3x3 matrix vectors representing the local x,y,z basis
    # and the 6x6 matrix of column vectors representing the local x,x',y,y',l,dp/p
    # phase space basis
    def BeamBasis(self):
        Ni = len(self.v)
        e3 = [self.v[0] / norm(self.v[0])]
        e2 = [np.cross(e3[0], np.array([0, -1, 0]))]
        e2[0] = e2[0] / norm(e2[0])
        e1 = [np.cross(e2[0], e3[0])]
        self.BasisM3 = [Basis3(e1[0], e2[0], e3[0])]
        self.BasisM6 = [Basis6(e1[0], e2[0], e3[0])]
        for i in range(1, Ni):
            e3.append(self.v[i] / norm(self.v[i]))
            e2.append(np.cross(e3[-1], e1[-1]))
            e2[-1] = e2[-1] / norm(e2[-1])
            e1.append(np.cross(e2[-1], e3[-1]))
            self.BasisM3.append(Basis3(e1[-1], e2[-1], e3[-1]))
            self.BasisM6.append(Basis6(e1[-1], e2[-1], e3[-1]))

    def PlotParticle(self, box=1.5, offsetX=0.5, offsetY=0.5, offsetZ=0):
        x = []
        y = []
        z = []
        R = []
        for i in range(len(self.r)):
            x.append(self.r[i][0])
            y.append(self.r[i][1])
            z.append(self.r[i][2])
            R.append(np.sqrt(x[-1]**2 + y[-1]**2))
        self.axs_3d.plot(x, y, z, color=self.LineColor,
                         linewidth=self.LineWidth)
        self.axs_3d.scatter(x[-1], y[-1], z[-1], s=15, c=self.LineColor)
        self.axs_pol.plot(R, z, color=self.LineColor,
                          linestyle=self.LineStyle, linewidth=self.LineWidth)
        self.axs_top.plot(x, y, color=self.LineColor,
                          linestyle=self.LineStyle, linewidth=self.LineWidth)
        self.axs_3d.set_xlim3d(-box / 2 + offsetX, box / 2 + offsetX)
        self.axs_3d.set_ylim3d(-box / 2 + offsetY, box / 2 + offsetY)
        self.axs_3d.set_zlim3d(-box / 2 + offsetZ, box / 2 + offsetZ)

# ------------------------------------------------------------------------------
# Plot 2D projection of trajectory
    def Plot2D(self, Type='poloidal'):
        x = []
        y = []
        z = []
        R = []
#		plt.figure(FIG)
        for i in range(len(self.r)):
            x.append(self.r[i][0])
            y.append(self.r[i][1])
            z.append(self.r[i][2])
            R.append(np.sqrt(x[-1]**2 + y[-1]**2))
        if Type == 'poloidal':
            PLOT = plt.plot(R, z, color=self.LineColor,
                            linestyle=self.LineStyle, linewidth=self.LineWidth)
        if Type == 'top':
            PLOT = plt.plot(x, y, color=self.LineColor,
                            linestyle=self.LineStyle, linewidth=self.LineWidth)
        return PLOT

# ------------------------------------------------------------------------------
# Initialize 3D Axes on figure
    def Figure3D(self, FIG=1):
        fig = plt.figure(FIG)
        ax = Axes3D(fig)
        return ax

# ------------------------------------------------------------------------------
# Plot trajectory in 3D
    def Plot3D(self, ax):
        x = []
        y = []
        z = []
        for i in range(len(self.r)):
            x.append(self.r[i][0])
            y.append(self.r[i][1])
            z.append(self.r[i][2])
        ax.plot(x, y, z, color=self.LineColor, linewidth=self.LineWidth)
        ax.scatter(x[-1], y[-1], z[-1], s=15, c=self.LineColor)
        return ax

    def Limits3D(self, ax, box=1.5, offsetX=0.5, offsetY=0.5, offsetZ=0):
        ax.set_xlim3d(-box / 2 + offsetX, box / 2 + offsetX)
        ax.set_ylim3d(-box / 2 + offsetY, box / 2 + offsetY)
        ax.set_zlim3d(-box / 2 + offsetZ, box / 2 + offsetZ)

    def FigureBV(self):
        self.fig_BV, self.axs_BV = plt.subplots(4, 3, figsize=(20, 8))
        self.axs_bx = self.axs_BV[0, 0]
        self.axs_bx.set_ylabel(r'Bx [T]')
        self.axs_bx.set_title('B-Field Components Along Trajectory')
        self.axs_by = self.axs_BV[1, 0]
        self.axs_by.set_ylabel(r'By [T]')
        self.axs_bz = self.axs_BV[2, 0]
        self.axs_bz.set_ylabel(r'Bz [T]')
        self.axs_b = self.axs_BV[3, 0]
        self.axs_b.set_ylabel(r'|B| [T]')
        self.axs_b.set_xlabel('S-coordinate [m]')

        self.axs_vx = self.axs_BV[0, 1]
        self.axs_vx.set_ylabel(r'$\beta_x$')
        self.axs_vx.set_title(
            r'Velocity Components Along Trajectory $\beta=v_i/c$')
        self.axs_vy = self.axs_BV[1, 1]
        self.axs_vy.set_ylabel(r'$\beta_y$')
        self.axs_vz = self.axs_BV[2, 1]
        self.axs_vz.set_ylabel(r'$\beta_z$')
        self.axs_v = self.axs_BV[3, 1]
        self.axs_v.set_ylabel(r'$\beta$')
        self.axs_v.set_xlabel('S-coordinate [m]')

        self.axs_bvx = self.axs_BV[0, 2]
        self.axs_bvx.set_ylabel(r'$\beta_x$')
        self.axs_bvy = self.axs_BV[1, 2]
        self.axs_bvy.set_ylabel(r'$\beta_y$')
        self.axs_bvz = self.axs_BV[2, 2]
        self.axs_bv = self.axs_BV[3, 2]
        self.axs_bv.set_ylabel(r'$\beta$')
        self.axs_bv.set_xlabel('S-coordinate [m]')
        self.axs_bv.set_ylabel(r'$\beta$')

    def PlotBV(self):
        Bx = []
        By = []
        Bz = []
        Bmag = []
        for i in range(len(self.B)):
            Bx.append(self.B[i][0])
            By.append(self.B[i][1])
            Bz.append(self.B[i][2])
            Bmag.append(norm(self.B[i]))
        self.axs_bx.plot(self.s, Bx)
        self.axs_by.plot(self.s, By)
        self.axs_bz.plot(self.s, Bz)
        self.axs_b.plot(self.s, Bmag)

        Vx = []
        Vy = []
        Vz = []
        Vmag = []
        for i in range(len(self.v)):
            vx, vy, vz = self.v[i] / cnt.c
            Vx.append(vx)
            Vy.append(vy)
            Vz.append(vz)
            Vmag.append(np.sqrt(vx**2 + vy**2 + vz**2))
        self.axs_vx.plot(self.s, Vx)
        self.axs_vy.plot(self.s, Vy)
        self.axs_vz.plot(self.s, Vz)
        self.axs_v.plot(self.s, Vmag)

        self.axs_bvx.plot(Vx, Bx)
        self.axs_bvy.plot(Vy, By)
        self.axs_bvz.plot(Vz, By)
        self.axs_bv.plot(Vmag, Bmag)

# ------------------------------------------------------------------------------
# Plot Magnetic Field components along beam trajectory
    def PlotB(self, FIG=2):
        Bx = []
        By = []
        Bz = []
        Bmag = []
        plt.figure(FIG)
        for i in range(len(self.B)):
            Bx.append(self.B[i][0])
            By.append(self.B[i][1])
            Bz.append(self.B[i][2])
            Bmag.append(norm(self.B[i]))
        plt.subplot(4, 1, 1)
        plt.plot(self.s, Bx)
        plt.ylabel(r'Bx [T]')
        plt.title('B-Field Components Along Trajectory')
        plt.subplot(4, 1, 2)
        plt.plot(self.s, By)
        plt.ylabel(r'By [T]')
        plt.subplot(4, 1, 3)
        plt.plot(self.s, Bz)
        plt.ylabel(r'Bz [T]')
        plt.subplot(4, 1, 4)
        plt.plot(self.s, Bmag)
        plt.ylabel(r'|B| [T]')
        plt.xlabel('S-coordinate [m]')

# ------------------------------------------------------------------------------
# Plot velocity components along beam trajectory
    def PlotV(self, FIG=3):
        Vx = []
        Vy = []
        Vz = []
        c0 = 2.998e8
        plt.figure(FIG)
        for i in range(len(self.v)):
            Vx.append(self.v[i][0] / c0)
            Vy.append(self.v[i][1] / c0)
            Vz.append(self.v[i][2] / c0)
        plt.subplot(3, 1, 1)
        plt.plot(self.s, Vx)
        plt.ylabel(r'$\beta_x$')
        plt.title(r'Velocity Components Along Trajectory $\beta=v_i/c$')
        plt.subplot(3, 1, 2)
        plt.plot(self.s, Vy)
        plt.ylabel(r'$\beta_y$')
        plt.subplot(3, 1, 3)
        plt.plot(self.s, Vz)
        plt.ylabel(r'$\beta_z$')
        plt.xlabel('S-coordinate [m]')

# ------------------------------------------------------------------------------
# Save magnetic field and curvature parameters
    def SaveFieldParameters(self, TFCurrent, Path='Output/'):
        # Save field and geometric parameters along trajectory
        np.savetxt(Path + 'Curvature_I_' +
                   str(int(TFCurrent)) + '.txt', self.k)
        np.savetxt(Path + 'SCoord_I_' + str(int(TFCurrent)) + '.txt', self.s)
        np.savetxt(Path + 'GradB_I_' +
                   str(int(TFCurrent)) + '.txt', self.gradB)
        np.savetxt(Path + 'GradBk_I_' +
                   str(int(TFCurrent)) + '.txt', self.gradBn)
        np.savetxt(Path + 'GradBn_I_' +
                   str(int(TFCurrent)) + '.txt', self.gradBk)


def Basis3(e1, e2, e3):
    Basis = np.matrix([
        [e1[0], e2[0], e3[0]],
        [e1[1], e2[1], e3[1]],
        [e1[2], e2[2], e3[2]]], float)
    return Basis


def Basis6(e1, e2, e3):
    Basis = np.matrix([
        [e1[0], 0.0, e2[0], 0.0, e3[0], 0.0],
        [0.0, e1[0], 0.0, e2[0], 0.0, e3[0]],
        [e1[1], 0.0, e2[1], 0.0, e3[1], 0.0],
        [0.0, e1[1], 0.0, e2[1], 0.0, e3[1]],
        [e1[2], 0.0, e2[2], 0.0, e3[2], 0.0],
        [0.0, e1[2], 0.0, e2[2], 0.0, e3[2]]], float)
    return Basis


# ==============================================================================
# ==============================================================================
# ======== Extra code for alternate trajectory integration methods =============
# ==============================================================================
# ==============================================================================

# Radius of Curvature method with perpendicular projection of B and constant dTheta = dS/R(B)
#BMag = 0.0
#vMag = 0.0
#hPara = np.zeros(3, float)
#hPerp = np.zeros(3, float)
# if False:
#    while (c1 or c2) and i < Nmax:
#
#        self.B.append(B.local(self.r[-1]) + Bv.local(self.r[-1]))
#
#        BMag = norm(self.B[-1])
#
#        vMag = norm(self.v[-1])
#
#        # parallel to velocity unit vector
#        hPara = self.v[-1] / vMag
#
#        # Vector along bending Radius
#        hRadius = np.cross(self.v[-1], self.B[-1])
#        hRadius = hRadius / norm(hRadius)
#
#        # perpendicular to B unit vector
#        hPerp = np.cross(hRadius, hPara)
#
#        # Magnitude of perpendicular projection of B
#        BPerp = np.dot(self.B[-1], hPerp)
#
#        # Cyclotron Frequency
#        Omega = self.q0 * BPerp / self.m0
#        dTheta = 0.001  # Omega*self.dt
#
#        # Larmor radius
#        rL = (self.m0 * vMag) / (self.q0 * BPerp)
#        print(rL, Omega, (rL * dTheta))
#        # Change in r
#
#        drPara = rL * np.sin(dTheta) * hPara
#
#        drRad = rL * (np.cos(dTheta) - 1.0) * hRadius
#
#        vPara = vMag * np.cos(dTheta) * hPara
#
#        vRad = vMag * np.sin(dTheta) * hRadius
#
#        self.r.append(self.r[-1] + drPara + drRad)
#
#        self.s.append(self.s[-1] + (rL * dTheta))
#
#        self.dS.append(self.s[-1] - self.s[-2])
#
#        self.a.append(qm * np.cross(self.v[-1], self.B[-1]))
#
#        self.v.append(vPara + vRad)
#
#        # Normalized Relativistic Parameters
#        self.Beta.append(self.v[-1] / c0)
#        self.beta.append(norm(self.Beta[-1]))
#        self.gamma.append(1.0 / (1.0 - self.beta[-1]**2))
#
#        # Check to see if beam np.crosses boundary
#        IN, NormalV, TangentV, IncidentV, RT, Xpol = Vessel.Xboundary(
#            self.r[-2], self.r[-1])
#
#        c1 = IN
#        c2 = self.s[-1] < Smin
#        i = i + 1
#        print(i)
##			self.Target = Target(NormalV,TangentV,IncidentV)
#
#    print('trajectory complete')
#    self.Target = Target(NormalV, TangentV, IncidentV,
#                         BFieldTF, BFieldVF, RT, Xpol)
#    self.BeamBasis()
#    print('Beam Coordinates Complete')
#    print(self.BasisM3)
# Boris Method with constant dTheta = dS/R(B)
# if False:
#    while (c1 or c2) and i < Nmax:
#
#        self.B.append(
#            np.array(B.local(self.r[-1])) + np.array(Bv.local(self.r[-1])))
#        BMag = norm(self.B[-1])
#
#        # parallel to B unit vector
#        hPara = self.B[-1] / BMag
#
#        # perpendicular to B unit vector
#        hPerp = (V - (V * hPara) * hPara)
#        hPerp = hPerp / norm(hPerp)
#
#        # Vector along bending Radius
#        hRadius = np.cross(hPara, hPerp)
#
#        # Cyclotron Frequency
#        Omega = self.q0 * BMag / self.m0
#        dTheta = Omega * self.dt
#
#        # Larmor radius
#        rL = (self.m0 * np.dot(self.v[-1], hPerp)) / (self.q0 * BMag)
#
#        # Change in r
#
#        drPara = self.dt * np.dot(self.v[-1], hPara) * hPara
#
#        drPerp = (rL * np.sin(dTheta)) * hPerp
#
#        drRad = rL * (np.cos(dTheta) - 1.0) * hRad
#
#        self.r.append(self.r[-1] + drPara + drPerp + drRad)
#
#        self.s.append(self.s[-1] + norm(self.r[-1] - self.r[-2]))
#
#        self.a.append(qm * np.cross(self.v[-1], self.B[-1]))
#
#        dv = np.dot(self.v[-1], hPara) * hPara + np.dot(self.v[-1],
#                                                        hPerp) * (np.cos(dTheta) * hPerp - np.sin(dTheta) * hRad)
#
#        self.v.append(self.v[-1] + dv)
#
#        # Normalized Relativistic Parameters
#        self.Beta.append(self.v[-1] / c0)
#        self.beta.append(norm(self.Beta[-1]))
#        self.gamma.append(1.0 / (1.0 - self.beta[-1]**2))
#
#        # Check to see if beam np.crosses boundary
#        IN, NormalV, TangentV, IncidentV, RT, Xpol = Vessel.Xboundary(
#            self.r[-2], self.r[-1])
#
#        c1 = IN
#        c2 = i * dS < Smin
#        i = i + 1
#        print(i)
##			self.Target = Target(NormalV,TangentV,IncidentV)
#
#    print('trajectory complete')
#    self.BeamBasis()
#    print('Beam Coordinates Complete')
#    self.Target = Target(NormalV, TangentV, IncidentV,
#                         BFieldTF, BFieldVF, RT, Xpol)
#    print('Target Complete')
#
#    self.NormalV = NormalV
#    self.IncidentV = IncidentV
#
