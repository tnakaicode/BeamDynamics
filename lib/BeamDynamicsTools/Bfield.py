import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sp
#from scipy.special import *

# Br = B0*R0/R * (1 + f(r,z,phi))
# Bz = BZ0
mu = 4 * np.pi * 1e-7

# ===============================================================================
# A collection B-Field Classes
# ===============================================================================

# ===============================================================================
# Simple 1/R toroidal B-Field Class
# ===============================================================================


class Bfield:
    # Br =  ra
    def __init__(self, B0, R0, B0z=0.0, fR=0, fz=0):
        self.B0 = B0
        self.B0z = B0z
        self.R0 = R0
        self.fR = fR
        self.fz = fz
        print('Bfield initialized')

    def localRZP(self, R, Z, Phi):
        Btor = self.B0 * self.R0 / R * (1 + self.fR)
        Bx = Btor * np.sin(Phi)
        By = Btor * np.cos(Phi)
        Bz = self.B0z
        return([Bx, By, Bz])

    def local(self, r):
        R = np.sqrt(r[0]**2 + r[1]**2)
        Phi = np.arctan(r[1] / r[0])
        Btor = self.B0 * self.R0 / R * (1 + self.fR)
        Bx = Btor * np.sin(Phi)
        By = Btor * np.cos(Phi)
        Bz = self.B0z
        if R > 1.3:
            B = np.array([0.0, 0.0, Bz])
        else:
            B = np.array([Bx, By, Bz])
        return B

# ===============================================================================
# Constant B-Field
# ===============================================================================


class Bfieldc:
    def __init__(self, B0=1.0, R0=1.0, B0z=0.0, fR=0, fz=0):
        self.B0 = B0
        self.B0z = B0z
        self.R0 = R0
        self.fR = fR
        self.fz = fz
        self.I0 = CalculateI0(self.B0, R0=0.66, NCoils=120.0)
        print('Bfield initialized')

    def local(self, r):
        Bx = 0.0
        By = self.B0
        Bz = 0.0
        return(np.array([Bx, By, Bz]))

# ===============================================================================
# Toroidal B-Field with discrete filament model
# ===============================================================================

# ======= Realistic TF Field ==================================================
# ======= Imported From BFieldDevelopment.py on 4/23/2013 =====================


class BfieldTF:
    # Generates Toroidal Field Coils
    def __init__(self, B0=1.0, R0=0.67, Phi0=2 * np.pi / 40, Ncoils=20, Rmin=0.1947965, Rmax=1.195229, Method='Filament'):

        TF = []
        self.B0 = B0
        self.R0 = R0
        self.Ncoils = Ncoils
        self.Rmin = Rmin
        self.Rmax = Rmax
        self.I0 = CalculateI0(self.B0, R0=0.67, NCoils=120.0)
        self.Method = Method

        for n in range(Ncoils):  # Outer Legs of TF
            TF.append(np.array([Rmax * np.cos(2 * np.pi * n / Ncoils + Phi0),
                                Rmax * np.sin(2 * np.pi * n / Ncoils + Phi0), -1.0]))

        # Inner Legs of TF -> np.array([ x , y , +/- direction ])
        for n in range(Ncoils):
            TF.append(np.array([Rmin * np.cos(2 * np.pi * n / Ncoils + Phi0),
                                Rmin * np.sin(2 * np.pi * n / Ncoils + Phi0), 1.0]))
        self.TF = TF

        # create np.array of TF coordinates to faster computation
        self.TFx = np.zeros(len(TF))
        self.TFy = np.zeros(len(TF))
        self.TFsign = np.zeros(len(TF))
        for i in range(len(TF)):
            self.TFx[i] = TF[i][0]
            self.TFy[i] = TF[i][1]
            self.TFsign[i] = TF[i][2]

    def PlotTF(self):
        plt.plt.figure(0)
        for n in range(len(TF)):
            plt.plot(TF[n][0], TF[n][1], 'ob')

    # Function that calculates Toroidal field at position R
    def local(self, RIN):
        if self.Method == 'Filament':
            Rx = RIN[0]
            Ry = RIN[1]
            AbsR = (Rx - self.TFx)**2 + (Ry - self.TFy)**2
            Bx = -1.0 * (self.B0 * self.R0 / self.Ncoils) * \
                (self.TFsign / AbsR) * (Ry - self.TFy)
            By = (self.B0 * self.R0 / self.Ncoils) * \
                (self.TFsign / AbsR) * (Rx - self.TFx)

            return np.array([sum(Bx), sum(By), 0.0])

        if self.Method == 'Filament0':
            R = np.array(RIN)
            B = np.array([0.0, 0.0, 0.0])
            Nc = len(self.TF) / 2  # /(2*np.pi)
            for n in range(len(self.TF)):
                AbsR = (R[0] - self.TF[n][0])**2 + (R[1] - self.TF[n][1])**2
                B[0] = B[0] - (self.B0 * self.R0 / Nc) * \
                    (self.TF[n][2] / AbsR) * (R[1] - self.TF[n][1])
                B[1] = B[1] + (self.B0 * self.R0 / Nc) * \
                    (self.TF[n][2] / AbsR) * (R[0] - self.TF[n][0])

#				if ( AbsR < 0.04**2 ):
#					B[0] = 0; B[1] = 0;
#					break
            return B

        if self.Method == 'Simple':
            R = np.sqrt(RIN[0]**2 + RIN[1]**2)
            B = np.array([0.0, 0.0, 0.0])
            r4 = 1.290259
            r3 = 1.100303  # coil Width [m]
            r2 = 0.302262
            r1 = 0.087331
            BTF = 0
            theta = np.arctan(RIN[1] / RIN[0])
            if R >= r1 and R <= r2:
                BTF = ((0.0 - (self.R0 * self.B0) / (r2)) / (r1 - r2)) * (R - r1)
            if R > r2 and R <= r3:
                BTF = (self.R0 * self.B0) / R
            if R > r3 and R <= r4:
                BTF = ((0.0 - (self.R0 * self.B0) / (r3)) / (r4 - r3)) * (R - r4)

            B = BTF * np.array([-np.sin(theta), np.cos(theta), 0.0])
            return B

# ===============================================================================
# Vertical B-Field generated with elliptic integrals
# ===============================================================================

# ======= Realistic VF Field ==================================================
# ======= Imported From BFieldDevelopment.py on 4/29/2013 =====================


class BfieldVF:

    def __init__(self, B0=None, I0=None, RCoil=[np.array([1.504188, 0.440817]), np.array([1.504188, -0.440817])]):

        self.RCoil = RCoil
        self.B0 = B0
        self.I0 = I0
        r0 = self.RCoil[0][0]
        self.r0 = r0
        self.nCoil = len(RCoil)

        if I0 == None and B0 != None:
            self.B0 = B0
            # /n0 #(2*r0)/(mu*self.B0)
            self.I0 = ((5.0 / 4.0)**(3.0 / 2.0) * (B0 * r0 / mu))

        if B0 == None and I0 != None:
            self.I0 = I0
            self.B0 = mu * I0 / (2.0 * r0)

    def local(self, R):
        #	RCoil=[np.array([1.0,0.0])]
        r = np.sqrt(R[0]**2 + R[1]**2)
        z1 = R[2]
        theta = np.arctan(R[1] / R[0])
        Br = 0.0
        Bz = 0.0
        BrR = 0.0
        for n in range(len(self.RCoil)):
            r0 = self.RCoil[n][0]
            z0 = self.RCoil[n][1]
            z = z1 - z0
            k = np.sqrt(4 * r * r0 / ((r + r0)**2 + z**2))
            IE = sp.ellipe(k)
            IK = sp.ellipk(k)

            Br = Br + (1.0 / r) * (z / np.sqrt((r + r0)**2 + z**2)) * \
                (-IK + IE * (r0**2 + r**2 + z**2) / ((r0 - r)**2 + z**2))

            Bz = Bz + (1.0 / np.sqrt((r + r0)**2 + z**2)) * \
                (IK + IE * (r0**2 - r**2 - z**2) / ((r0 - r)**2 + z**2))

            if ((r - r0)**2 + z**2 < 0.1**2):
                Br = 0
                Bz = 0
                break

        BV = (mu * self.I0) / (2.0 * r0) / self.nCoil * \
            np.array([Br * np.cos(theta), Br * np.sin(theta), Bz])
        return BV


# =============================================================================
# ======= Function to convert I into B0 =======================================
# =============================================================================

def CalculateB0(I0, R0=0.66, NCoils=120.0):
    mu = 4.0e-7 * np.pi
    B0 = (mu / (2 * np.pi)) * I0 * NCoils / R0
    return B0


def CalculateI0(B0, R0=0.66, NCoils=120.0):
    mu = 4.0e-7 * np.pi
    I0 = ((2 * np.pi) / mu) * (R0 / NCoils) * B0
    return I0
