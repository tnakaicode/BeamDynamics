import pylab as pl
import matplotlib.pyplot as plt
from pylab import det
import numpy as np
import math
from numpy.linalg import inv


class Ellipse:
    def __init__(self, SIG):
        self.Sigma = SIG
# ------------------------------------------------------------------------------
# 2x2 sigma matrices
        self.SigX = self.Sigma[0:2, 0:2]
        self.SigY = self.Sigma[2:4, 2:4]
        self.SigZ = self.Sigma[4:6, 4:6]

# ------------------------------------------------------------------------------
# Emittance in each phase place
        # self.EpsilonX = self.Sigma[0,0]*self.Sigma[1,1] - self.Sigma[0,1]**2
        self.EmittenceX = np.sqrt(math.fabs(det(self.SigX)))
        self.EmittenceY = np.sqrt(math.fabs(det(self.SigY)))
        self.EmittenceZ = np.sqrt(math.fabs(det(self.SigZ)))

# ------------------------------------------------------------------------------
# Ellipse (Twiss) parameters in each phase plane
        self.TwissXX1 = np.array([-(self.SigX[0, 1]), self.SigX[0, 0],
                               self.SigX[1, 1], self.EmittenceX**2] / self.EmittenceX)
        self.TwissYY1 = np.array([-(self.SigY[0, 1]), self.SigY[0, 0],
                               self.SigY[1, 1], self.EmittenceY**2] / self.EmittenceY)
        self.TwissZZ1 = np.array([-(self.SigZ[0, 1]), self.SigZ[0, 0],
                               self.SigZ[1, 1], self.EmittenceZ**2] / self.EmittenceZ)

# ------------------------------------------------------------------------------
# Beam's spatial width in each phase plane
        self.WidthX = np.sqrt(self.TwissXX1[1] * self.TwissXX1[3])
        self.WidthY = np.sqrt(self.TwissYY1[1] * self.TwissYY1[3])
        self.WidthZ = np.sqrt(self.TwissZZ1[1] * self.TwissZZ1[3])

# ------------------------------------------------------------------------------
# Beam's angular envelope in each phase plane
        self.DivergenceX = np.sqrt(self.TwissXX1[2] * self.TwissXX1[3])
        self.DivergenceY = np.sqrt(self.TwissYY1[2] * self.TwissYY1[3])
        self.DivergenceZ = np.sqrt(self.TwissZZ1[2] * self.TwissZZ1[3])

# ------------------------------------------------------------------------------
# Ellipse parameter in the transverse spatial plane
        self.EmittenceXY = np.sqrt(
            det(np.matrix([[SIG[0, 0], SIG[0, 2]], [SIG[2, 0], SIG[2, 2]]])))
        self.TwissXY = np.array(
            [-SIG[0, 2], SIG[0, 0], SIG[2, 2], self.EmittenceXY**2]) / self.EmittenceXY

        self.EmittenceXZ = np.sqrt(
            det(np.matrix([[SIG[0, 0], SIG[0, 4]], [SIG[4, 0], SIG[4, 4]]])))
        self.TwissXZ = np.array(
            [-SIG[0, 4], SIG[0, 0], SIG[4, 4], self.EmittenceXZ**2]) / self.EmittenceXZ

        self.EmittenceYZ = np.sqrt(
            det(np.matrix([[SIG[2, 2], SIG[2, 4]], [SIG[4, 2], SIG[4, 4]]])))
        self.TwissYZ = np.array(
            [-SIG[2, 4], SIG[2, 2], SIG[4, 4], self.EmittenceYZ**2]) / self.EmittenceYZ

# ------------------------------------------------------------------------------
# Plotting Attributes
        self.LineColor = 'r'
        self.LineWidth = 1
        self.LineStyle = '-'

# ------------------------------------------------------------------------------
# Return Spatial Width
    def SpatialWidth(self):
        return self.WidthX, self.WidthY, self.WidthZ

# ------------------------------------------------------------------------------
# Return Angular Width
    def AngularWidth(self):
        return self.DivergenceX, self.DivergenceY, self.DivergenceZ

# ------------------------------------------------------------------------------
# Generate points along an ellipse a give set of twiss parameters
    def GenerateXY(self, TWISS, NPoints=1000):
        Theta = np.linspace(0, 2 * np.pi, NPoints)
        XPoints = np.zeros((NPoints), float)
        YPoints = np.zeros((NPoints), float)
        m11 = np.sqrt(math.fabs(TWISS[1]))
        m21 = -TWISS[0] / np.sqrt(math.fabs(TWISS[1]))
        m22 = 1 / np.sqrt(math.fabs(TWISS[1]))
        Radius = np.sqrt(math.fabs(TWISS[3]))
        m12 = 0
        PHI = np.arctan(2.0 * TWISS[0] / (TWISS[2] - TWISS[1])) / 2.0
        for i in range(NPoints):
            XPoints[i] = Radius * (m11 * np.cos(Theta[i]) + m12 * np.sin(Theta[i]))
            YPoints[i] = Radius * (m21 * np.cos(Theta[i]) + m22 * np.sin(Theta[i]))
        return XPoints, YPoints

# ------------------------------------------------------------------------------
# Calculate Mismatch factor between self and another ellipse E1
    def MismatchFactor(self, E1, Type=1):

        def MFactor(Twiss0, Twiss1, Type):
            R = Twiss0[1] * Twiss1[2] + Twiss0[2] * \
                Twiss1[1] - 2.0 * Twiss0[0] * Twiss1[0]
            M = (0.5 * (R + np.sqrt(R**2 - 4.0)))**(0.5 * Type) - 1.0
            return M

        Mx = MFactor(self.TwissXX1, E1.TwissXX1, Type)
        My = MFactor(self.TwissYY1, E1.TwissYY1, Type)
        Mz = MFactor(self.TwissZZ1, E1.TwissZZ1, Type)
        Mxy = MFactor(self.TwissXY, E1.TwissXY, Type)

        return [Mx, My, Mz, Mxy]


#	def ProjectXY(XPoints,YPoints,Axz,Ayz): #Axz = angle in XZ plane
#		for i in range(len(XPoints)):
#			XPoints[i] = XPoints[i]/np.cos(Axz)/np.sin(Ayz)
#			YPoints[i] = YPoints[i]/np.cos(Ayz)/np.sin(Axz)
#		return XPoints,YPoints
#
#	def PlotProjectionXY(self,Axz=0,Ayz=0,FIG=0,NPoints=1000,Mod='-',Title = ' ',Label=''):
# f=plt.figure(FIG)
#		X,Y = self.GenerateXY(self.TwissXY,NPoints)
#		X = X/np.cos(Axz)
#		Y = Y/np.cos(Ayz)
#		plt.plot(X,Y,Mod,label=Label); plt.xlabel('X [mm]');  plt.ylabel('Y [mm]');

# ------------------------------------------------------------------------------
# Plot transverse spatial projection

    def PlotXY(self, NPoints=1000, L=30.0, Mod='-', Label='', Title=' ', Scale=1.0, Rotate=False):
        X, Y = self.GenerateXY(self.TwissXY, NPoints)
        if Rotate == True:
            Y, X = self.GenerateXY(self.TwissXY, NPoints)
            Y = (-1) * Y
        PLOT = plt.plot(Scale * X, Scale * Y, Mod, label=Label)
        plt.xlabel('X [mm]')
        plt.ylabel('Y [mm]')
        L = max([max(X), max(Y)]) * 1.2
        plt.xlim(-L, L)
        plt.ylim(-L, L)
#		plt.xlim(-L,L); plt.ylim(-L,L)
        return PLOT

# ------------------------------------------------------------------------------
# Plot X,X' phase plane projection
    def PlotXX1(self, NPoints=1000, L=30.0, Mod='-', Label='', Title=' ', Scale=1.0):
        X, X1 = self.GenerateXY(self.TwissXX1, NPoints)
        PLOT = plt.plot(Scale * X, Scale * X1, Mod, label=Label)
        plt.xlabel('X [mm]')
        plt.ylabel(r'X$^\prime$ [mrad]')
        L = max([max(X), max(X1)]) * 1.2
        plt.xlim(-L, L)
        plt.ylim(-L, L)
        return PLOT

# ------------------------------------------------------------------------------
# Plot Y,Y' phase plane projection
    def PlotYY1(self, NPoints=1000, L=30.0, Mod='-', Label='', Title=' ', Scale=1.0):
        Y, Y1 = self.GenerateXY(self.TwissYY1, NPoints)
        PLOT = plt.plot(Scale * Y, Scale * Y1, Mod, label=Label)
        plt.xlabel('Y [mm]')
        plt.ylabel(r'Y$^\prime$ [mrad]')
        L = max([max(Y), max(Y1)]) * 1.2
        plt.xlim(-L, L)
        plt.ylim(-L, L)
        return PLOT

# ------------------------------------------------------------------------------
# Plot Z,Z' phase plane projection
    def PlotZZ1(self, NPoints=1000, L=30.0, Mod='-', Label='', Title=' ', Scale=1.0):
        Z, Z1 = self.GenerateXY(self.TwissZZ1, NPoints)
        PLOT = plt.plot(Scale * Z, Scale * Z1, Mod, label=Label)
        plt.xlabel(r'$\ell$ [mm]')
        plt.ylabel(r'$\Delta$P/P [mrad]')
        L = max([max(Z), max(Z1)]) * 1.2
        plt.xlim(-L, L)
        plt.ylim(-L, L)
        return PLOT

# ------------------------------------------------------------------------------
# Project transverse beam spot onto off-normal surface
    def ProjectOffNormal(self, SigmaBasis, TargetBasis, Scale=1.0, Label='', Title=' ', NPoints=1000, Mod='-'):
        X, Y = self.GenerateXY(self.TwissXY, NPoints)
        Bs = np.matrix(SigmaBasis[:, [0, 2]])
        Bt = np.matrix(TargetBasis[:, [0, 2]])
        Ms = Scale * np.eye(2)
        Xp = []
        Yp = []
        Bdot = Bt * Bs.T
        Bproj = np.matrix(np.zeros((2, 2), float))
        for i in [0, 1]:
            for j in [0, 1]:
                Bproj[i, j] = 1.0 / (Bdot[i, j])
        b11 = 1.0 / (Bt[:, 0].T * Bs[:, 0])**2
        b22 = 1.0 / (Bt[:, 1].T * Bs[:, 1])**2
        Bproj = np.matrix([[b11[0, 0], 0.0], [0.0, b22[0, 0]]])
#		Ang =
#		print Bproj
        for i in range(len(X)):
            V = np.transpose(np.matrix([X[i], Y[i]]))
#			Vp = Bproj * (Bt.T * Bs) * Bs * V
#			Vp = Mrot * Bproj *V
#			Vp =Bt.T * Bs * V
            Vp = Ms * inv(Bs.T * Bt) * V
            Xp.append(Vp[0, 0])
            Yp.append(-Vp[1, 0])
        self.ProjectionX = np.array(Xp) / 1000.0
        self.ProjectionY = np.array(Yp) / 1000.0
        return self.ProjectionX, self.ProjectionY

    def PlotProjectionXY(self, offsetX=0.0, offsetY=0.0, Mod='-', Label=''):
        Xp = self.ProjectionX + offsetX
        Yp = self.ProjectionY + offsetY
        plt.plot(Xp, Yp, Mod, label=Label, color=self.LineColor)
        plt.xlabel('X [mm]')
        plt.ylabel('Y [mm]')

    def PlotProjectionPolPhi(self, Mod='-', Label=''):  # NOT FINISHED
        Xp = self.ProjectionX
        Yp = self.ProjectionY
        plt.plot(Xp, Yp, Mod, label=Label)
        plt.xlabel('Poloidal Position [mm]')
        plt.ylabel('Toroidal Angle [degrees]')


# ------------------------------------------------------------------------------
# Save XY projection points

    def PrintProjection(self, FileName='ProjectionXY'):
        Output = np.transpose(np.vstack((self.ProjectionX, self.ProjectionY)))
        np.savetxt(FileName, Output)

# ------------------------------------------------------------------------------
# Plot All projections
    def PlotALL(self, FIG=0, NPoints=1000, Mod='-', Title=' '):

        f = plt.figure(FIG)
        f.text(.5, .95, Title, horizontalalignment='center')

        X, Y = self.GenerateXY(self.TwissXX1, NPoints)
        plt.subplot(2, 3, 1)
        plt.plot(X, Y, Mod)
        plt.xlabel('X [mm]')
        plt.ylabel(r'$\Delta$Px / P [mrad]')

        X, Y = self.GenerateXY(self.TwissYY1, NPoints)
        plt.subplot(2, 3, 2)
        plt.plot(X, Y, Mod)
        plt.xlabel('Y [mm]')
        plt.ylabel(r'$\Delta$Py / P [mrad]')

        X, Y = self.GenerateXY(self.TwissZZ1, NPoints)
        plt.subplot(2, 3, 3)
        plt.plot(X, Y, Mod)
        plt.xlabel('Z [mm]')
        plt.ylabel(r'$\Delta$Pz / P [mrad]')

#		plt.figure(FIG+1)
        X, Y = self.GenerateXY(self.TwissXY, NPoints)
        plt.subplot(2, 3, 4)
        plt.plot(X, Y, Mod)
        plt.xlabel('X [mm]')
        plt.ylabel('Y [mm]')
        L = 20
        plt.xlim(-L, L)
        plt.ylim(-L, L)

        X, Y = self.GenerateXY(self.TwissXZ, NPoints)
        plt.subplot(2, 3, 5)
        plt.plot(X, Y, Mod)
        plt.xlabel('X [mm]')
        plt.ylabel('Z [mm]')
        L = 20
        plt.xlim(-L, L)
        plt.ylim(-L, L)

        X, Y = self.GenerateXY(self.TwissYZ, NPoints)
        plt.subplot(2, 3, 6)
        plt.plot(X, Y, Mod)
        plt.xlabel('Y [mm]')
        plt.ylabel('Z [mm]')
        L = 20
        plt.xlim(-L, L)
        plt.ylim(-L, L)

        plt.subplots_adjust(hspace=0.35)
        plt.subplots_adjust(wspace=0.35)

#		plt.xlim([-2,2]); plt.ylim([-2,2])
