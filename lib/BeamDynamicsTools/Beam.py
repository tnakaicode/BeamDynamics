# Beam.py

import numpy as np
#import scipy as sp
import pylab as pl
from numpy.linalg import inv, norm
from .Trajectory import *
from .Target import *
from .Ellipse import *


class Beam(Trajectory):
    # inputs:
    # sigma = 6x6 sigma matrix
    # s0 = 3x3 matrix for local beam coordinate system
    def __init__(self, trajectory, sigma0):
        self.sigma = [sigma0]

        self.q0 = trajectory.q0
        self.c0 = trajectory.c0
        self.m0 = trajectory.m0
        self.I0 = trajectory.I0
        self.Frequency = trajectory.Frequency
#		self.Z0 = trajectory.z0
        self.A0 = trajectory.A0
#		self.Beta = trajectory.Beta
#		self.Ellipse = [Ellipse(sigma0)]
        self.r = trajectory.r
        self.v0 = trajectory.v0
        self.v = trajectory.v
        self.Beta = trajectory.Beta
        self.beta = trajectory.beta
        self.gamma = trajectory.gamma
        self.a = trajectory.a
        self.B = trajectory.B
        self.s = trajectory.s
        self.k = trajectory.k
        self.Rc = trajectory.Rc
        self.dS = trajectory.dS
        self.dt = trajectory.dt
        self.BasisM3 = trajectory.BasisM3
        self.BasisM6 = trajectory.BasisM6
        self.gradB = trajectory.gradB
        self.gradBk = trajectory.gradBk
        self.gradBn = trajectory.gradBn
        self.gradBx = trajectory.gradBx
        self.gradBy = trajectory.gradBy
        self.NormalV = trajectory.target.NormalV
        self.IncidentV = trajectory.target.IncidentV
        self.target = trajectory.target

        # lists of matrices for reverse calculations
        self.Revsigma = []
        self.RevTransferM = []

#		self.e1=trajectory.e1
#		self.e2=trajectory.e2
#		self.e3=trajectory.e3

# ===============================================================================
# Calculate evolution of sigma np.matrix along trajectory
# ===============================================================================

    def Trace(self):
        Ni = len(self.r)
        self.TransferM = []
        for i in range(Ni):
            #			D = self.Drift(self.dS[i])
            # Mb is the matrix form of Acc = (q/m) v x B
            S = self.BasisM3[i]  # np.matrix(identity(6))
#			print self.B[i]
#			B = self.BMatrix(S,self.B[i],self.dS[i])
            M = self.BMatrix(i)
#			Q = self.SpaceCharge(i,self.sigma[i]) #,self.dS[i],self.beta[i],self.gamma[i])
#			M = Q*M
            self.TransferM.append(M)
            self.sigma.append(M * self.sigma[-1] * M.T)
            self.sigma[-1] = self.SpaceCharge(i, self.sigma[-1])
            # print i

        self.target.sigma = self.sigma[-1]
        self.target.BeamBasis = self.BasisM6[-1]
        self.target.ProjEllipse = Ellipse(self.sigma[-1])
        self.target.ProjectionX, self.target.ProjectionY = self.target.ProjEllipse.ProjectOffNormal(
            self.target.BeamBasis, self.target.TargetBasis)
#		self.target.ProjectSigma(self,SB=self.target.BeamBasis,ST=self.TargetBasis)


# ===============================================================================
# Calculate reverse evolution of sigma matrix starting at target
# ===============================================================================

    def ReverseTrace(self, sigmaF):
        Ni = len(self.r)
        RevTransferM = []
        Revsigma = [sigmaF]
        for i in range(Ni - 1, -1, -1):
            D = self.Drift(self.dS[i])
            # Mb is the matrix form of Acc = (q/m) v x B
            S = self.BasisM3[i]  # np.matrix(identity(6))
            B = self.BMatrix(S, self.B[i], self.dS[i])
            M = inv(B)
            RevTransferM.append(M)
            Revsigma.append(M * Revsigma[-1] * M.T)
            print(i)
        self.Revsigma.append(Revsigma)
        self.RevTransferM.append(RevTransferM)

    def Drift(self, ds=1e-3):
        Mdrift = np.matrix([
            [1, ds, 0, 0, 0, 0],
            [0, 1, 0, 0, 0, 0],
            [0, 0, 1, ds, 0, 0],
            [0, 0, 0, 1, 0, 0],
            [0, 0, 0, 0, 1, ds],
            [0, 0, 0, 0, 0, 1]], float)
        # print Mdrift
        return Mdrift


# ===============================================================================
# Generates transfer matrices for magnetic steering
# ===============================================================================

    def BMatrix(self, IND):
        #		print 'Index = '+str(IND)
        #		Bperp = matrix([[np.dot(self.B[IND],BasisM3[IND][:,0])],[np.dot(self.B[IND],BasisM3[IND][:,1])],[0.0]])
        # ------------------------------------------------------------------------------
        # project B onto local coordinates
        Bx = np.dot(self.B[IND], self.BasisM3[IND][:, 0])
        By = np.dot(self.B[IND], self.BasisM3[IND][:, 1])
        dS = self.dS[IND]
        G2 = self.gamma[IND]**2

        QP = self.q0 / (self.m0 * self.v0)
# ------------------------------------------------------------------------------
# Calculate field gradients in local coordinates

#		kx = 0.0*QP * (self.gradBx[IND]);
        kx = QP * (self.gradBx[IND])
        Kappax = QP * Bx
        Kx = (kx + Kappax**2)  # /1000.0;
        Ax = np.sqrt(abs(Kx))

#		ky = QP * (self.gradBx[IND]);
        ky = QP * (self.gradBy[IND])
        Kappay = QP * By
        Ky = (ky + Kappay**2)
        # print kx,ky#'Ky %0.000f' %Ky   #/1000.0; print Ky
        Ay = np.sqrt(abs(Ky))

# ------------------------------------------------------------------------------
# Generate matrix entries based on linear model
        if norm(self.B[IND]) != 0:
            if Kx == 0:
                Cx = 1.0
                Cx1 = 0.0
                Sx = dS
                Sx1 = 1.0
                Dx = 0.0
                Dx1 = 0.0

            if Kx > 0:
                Cx = np.cos(Ax * dS)
                Cx1 = - Ax * np.sin(Ax * dS)
                Sx = (1.0 / Ax) * np.sin(Ax * dS)
                Sx1 = np.cos(Ax * dS)
                Dx = (1.0 / Ax) * (1.0 - np.cos(Ax * dS))
                Dx1 = np.sin(Ax * dS)

            if Kx < 0:
                Cx = np.cosh(Ax * dS)
                Cx1 = - Ax * np.sinh(Ax * dS)
                Sx = (1.0 / Ax) * np.sinh(Ax * dS)
                Sx1 = np.cosh(Ax * dS)
                Dx = (1.0 / Ax) * (1.0 - np.cosh(Ax * dS))
                Dx1 = np.sinh(Ax * dS)

            if Ky == 0:
                Cy = 1.0
                Cy1 = 0.0
                Sy = dS
                Sy1 = 1.0
                Dy = 0.0
                Dy1 = 0.0

            if Ky > 0:
                Cy = np.cos(Ay * dS)
                Cy1 = - Ay * np.sin(Ay * dS)
                Sy = (1.0 / Ay) * np.sin(Ay * dS)
                Sy1 = np.cos(Ay * dS)
                Dy = (1.0 / Ay) * (1.0 - np.cos(Ay * dS))
                Dy1 = np.sin(Ay * dS)

            if Ky < 0:
                Cy = np.cosh(Ay * dS)
                Cy1 = - Ay * np.sinh(Ay * dS)
                Sy = (1.0 / Ay) * np.sinh(Ay * dS)
                Sy1 = np.cosh(Ay * dS)
                Dy = (1.0 / Ay) * (1.0 - np.cosh(Ay * dS))
                Dy1 = np.sinh(Ay * dS)

# ------------------------------------------------------------------------------
# Populate transfer matrices with calculated entry values
            Mb = np.matrix([
                [Cx, Sx, 0, 0, 0, Dx],
                [Cx1, Sx1, 0, 0, 0, Dx1],
                [0, 0, Cy, Sy, 0, Dy],
                [0, 0, Cy1, Sy1, 0, Dy1],
                [-Dx1, -Dx, -Dy1, -Dy, 1, dS / G2 / 1.],
                [0, 0, 0, 0, 0, 1]], float)

# ------------------------------------------------------------------------------
# If there is no field, use tranfer matrix for field free drift
        else:
            Mb = np.matrix([
                [1, dS, 0, 0, 0, 0],
                [0, 1, 0, 0, 0, 0],
                [0, 0, 1, dS, 0, 0],
                [0, 0, 0, 1, 0, 0],
                [0, 0, 0, 0, 1, dS / G2],
                [0, 0, 0, 0, 0, 1]], float)

        return Mb

# =============================================================================
# =========== Space Charge Effects ============================================
# =============================================================================
    def SpaceCharge(self, IND, sigmaIN):
        G2 = self.gamma[IND]**2

        sigma = sigmaIN

# ------------------------------------------------------------------------------
# Rotate upright in XY Plane
        ThetaXY = 0.5 * np.arctan(2 * sigma[0, 2] / (sigma[2, 2] - sigma[0, 0]))
        if np.isnan(ThetaXY):
            ThetaXY = 0.0
        C = np.cos(ThetaXY)
        S = np.sin(ThetaXY)
        Rxy = np.matrix([
            [C, 0, -S, 0, 0, 0],
            [0, C, 0, -S, 0, 0],
            [S, 0, C, 0, 0, 0],
            [0, S, 0, C, 0, 0],
            [0, 0, 0, 0, 1, 0],
            [0, 0, 0, 0, 0, 1]], float)
        sigma = Rxy * sigma * Rxy.T

# ------------------------------------------------------------------------------
# Rotate upright in YZ Plane (Most Important)
        ThetaYZ = -0.5 * np.arctan(2 * sigma[2, 4] / (sigma[4, 4] - sigma[2, 2]))
        if np.isnan(ThetaYZ):
            ThetaYZ = 0.0
        C = np.cos(ThetaYZ)
        S = np.sin(ThetaYZ)
        Ryz = np.matrix([
            [1, 0, 0, 0, 0, 0],
            [0, 1, 0, 0, 0, 0],
            [0, 0, C, 0, -S, 0],
            [0, 0, 0, C, 0, -S],
            [0, 0, S, 0, C, 0],
            [0, 0, 0, S, 0, C]], float)
        sigma = Ryz * sigma * Ryz.T

# ------------------------------------------------------------------------------
# Rotate upright in XZ Plane
        C = 1.0
        S = 0.0
        ThetaZX = -0.5 * np.arctan(2 * sigma[4, 0] / (sigma[0, 0] - sigma[4, 4]))
        if np.isnan(ThetaZX):
            ThetaZX = 0.0
        C = np.cos(ThetaZX)
        S = np.sin(ThetaZX)
        Rzx = np.matrix([
            [C, 0, 0, 0, -S, 0],
            [0, C, 0, 0, 0, -S],
            [0, 0, 1, 0, 0, 0],
            [0, 0, 0, 1, 0, 0],
            [S, 0, 0, 0, C, 0],
            [0, S, 0, 0, 0, C]], float)
        sigma = Rzx * sigma * Rzx.T

# ------------------------------------------------------------------------------
# Rotation Matrix
        Rotate = Rxy * Ryz * Rzx
#		sigma = Rotate * sigma * Rotate.T

        # print sigma[0,2],sigma[2,4], sigma[4,0]
#		print ThetaXY, ThetaYZ, ThetaZX
#		print 0.5 * np.arctan(2*sigma[0,2]/(sigma[2,2]-sigma[0,0])), 0.5 * np.arctan(2*sigma[2,4]/(sigma[4,4]-sigma[2,2])),-0.5 * np.arctan(2*sigma[4,0]/(sigma[0,0]-sigma[4,4]))
#		print sigma[0,2],sigma[2,4],sigma[4,0]
        # Beam semiaxes
        Const = 1.0e-3 * np.sqrt(5.0)

        rx = np.sqrt(abs(sigma[0, 0])) * Const
        ry = np.sqrt(abs(sigma[2, 2])) * Const
        rz = np.sqrt(abs(sigma[4, 4])) * Const / self.gamma[IND]

#		sigma = (1.0/5.0)*sigma

#		print rx,ry,rz

# ------------------------------------------------------------------------------
# calculate normalized radial component for form factor fit
        p = self.gamma[IND] * rz / np.sqrt(rx * ry)

        # Form factor f polyfit Coefficient
        f = 0.0
        if p < 1.0:
            C0 = [0.32685993, -1.10422029, 1.64157723, -1.52987752, 0.99919503]
            f = np.polyval(C0, p)
        if p >= 1.0:
            C0 = [0.39122685, -0.97242357, 0.74316934, 0.1744343, -0.0020169]
# ------------------------------------------------------------------------------
# Calculate space charge form factor f from polynomial fit
            f = np.polyval(C0, 1.0 / p)

# ------------------------------------------------------------------------------
# Calculate Constants for E-Field Calculation
        k = 1.0 / (4 * np.pi * 8.85e-12)
        Lambda = self.c0 / self.Frequency
        Q = 3.0 * self.I0 * Lambda / self.c0

# ------------------------------------------------------------------------------
# Calculate Space charge E-Field Components
        Ex = (k * Q / self.gamma[IND]**2) * (1.0 - f) / (rx * (rx + ry) * rz)
        Ey = (k * Q / self.gamma[IND]**2) * (1.0 - f) / (ry * (rx + ry) * rz)
        Ez = (k * Q) * f / (rx * ry * rz)  # /np.sqrt(5)

# ------------------------------------------------------------------------------
# Constant to convert E-field to delta Xi'
        d = (self.q0 * self.dS[IND]) / (self.m0 * self.c0 **
                                        2 * self.beta[IND]) * 1e3 * np.sqrt(1 / 5.0)
#		print d*Ez

# ------------------------------------------------------------------------------
# Apply SPace Charge Impulses to momentum
        ME = np.matrix([
            [1, 0, 0, 0, 0, 0],
            [d * Ex, 1, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0],
            [0, 0, d * Ey, 1, 0, 0],
            [0, 0, 0, 0, 1, 0],
            [0, 0, 0, 0, d * Ez, 1]], float)

#		print d*Ex/np.sqrt(sigma[1,1]),d*Ey/np.sqrt(sigma[3,3]),d*Ez/np.sqrt(sigma[5,5])
        # Rotate back to orignial orientation
##		sigma = (Rxy.T * Ryz.T * Rzx.T) * sigma * (Rzx * Ryz * Rxy)
#		return (Rzx.T * Ryz.T * Rxy.T) * ME * sigma * ME.T * ( Rxy *Ryz * Rzx )
#		return (Rzx * Ryz * Rxy) * ME * sigma * ME.T * (Rzx.T * Ryz.T * Rxy.T)
#		return ME * sigma * ME.T
        return Rotate.T * ME * sigma * ME.T * Rotate

# =============================================================================
# =============================================================================
# =============================================================================
# =========== Old Code Below ==================================================
# =============================================================================
# =============================================================================
# =============================================================================

    def BMatrix1(self, Basis, Bin, dS):

        Bperp = np.matrix([[np.dot(Bin, Basis[:, 0])], [
                       np.dot(Bin, Basis[:, 1])], [0.0]])

        B0 = norm(Bperp)
        if B0 != 0:
            r = (self.m0 * self.v0) / (self.q0 * B0)
            da = dS / r
            C = np.cos(da)
            S = np.sin(da)

            dA = np.arctan(Bperp[0, 0] / Bperp[1, 0])
            Cr = np.cos(dA)
            Sr = np.sin(dA)

            R0 = np.matrix([
                [Cr, 0, Sr, 0, 0, 0],
                [0, Cr, 0, Sr, 0, 0],
                [-Sr, 0, Cr, 0, 0, 0],
                [0, -Sr, 0, Cr, 0, 0],
                [0, 0, 0, 0, 1, 0],
                [0, 0, 0, 0, 0, 1]], float)

            Mb = np.matrix([
                [1, dS, 0, 0, 0, 0],
                [0, 1, 0, 0, 0, 0],
                [0, 0, C, r * S, 0, r - r * C],
                [0, 0, -S / r, C, 0, S],
                [0, 0, 0, 0, 1, dS],
                [0, 0, 0, 0, 0, 1]], float)

            Mb = R0.T * Mb * R0
#			Mb = R0 * Mb * R0.T
        else:
            Mb = self.Drift(self.dS[IND])

        return Mb

    def BMatrix0(self, Vin, Bin):
        Bn = (self.q0 / self.m0 / self.v0) * Bin  # / 1e3 #(self.dt/self.v0)
        Vn = Vin / self.v0 * self.dS

        Fxy = Vn[0] * Bn[1]
        Fyx = -Vn[1] * Bn[0]
        Fyz = Vn[1] * Bn[2]
        Fzy = -Vn[2] * Bn[1]
        Fzx = Vn[2] * Bn[0]
        Fxz = -Vn[0] * Bn[2]

        Mb = np.matrix([
            [1, 0, 0, 0, 0, 0],
            [0, 1, 0, Fyz, 0, Fyx],
            [0, 0, 1, 0, 0, 0],
            [0, Fxz, 0, 1, 0, Fzx],
            [0, 0, 0, 0, 1, 0],
            [0, Fxy, 0, Fyx, 0, 1]], float)
#		Mb = matrix(identity(6))
        return Mb

#		sigma = ME * sigma * ME.T
#		print d*Ex,d*Ey,d*Ez

#		print Rxy
#		print Ryz
#		print Rzx
#		print ME
#		print sigma
#		print ME
#		print ThetaXY, ThetaYZ, ThetaZX
