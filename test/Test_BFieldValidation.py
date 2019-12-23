import sys
sys.path.append('../lib/BeamDynamicsTools/')
from Bfield import *

import numpy as np
import pylab as pl
from scipy.special import *
from numpy.linalg import norm
import matplotlib as mpl

R = np.array([1, 0, 0])
RInj = [1.798, -0.052, 0.243]

# ------------------------------------------------------------------------------
# Define Toroidal Field
TF = BfieldTF()

# ------------------------------------------------------------------------------
# Define Vertical Field
VF = BfieldVF(B0=1.0)

# ===============================================================================
# B Field Plotting Functions
# ===============================================================================

# ------------------------------------------------------------------------------
# Plot Elliptic Integrals of the 1st and 2nd kind, used for vertical field


def PlotEllipticInt():
    Ni = 1000
    Ek = []
    Ee = []
    for i in range(Ni):
        K1 = ellipk(1.0 * i / Ni)
        print(1.0 * i / Ni)  # K1
        Ek.append(K1)  # print Ek[-1]
        E1 = ellipe(1.0 * i / Ni)
        Ee.append(E1)  # print Ee[-1]
    plt.plot(Ek)
    plt.plot(Ee)

# ------------------------------------------------------------------------------
# Plot Toroidal field in 1D


def PlotTF1D(TF):
    Rmin = 0.1947965
    Rmax = 1.195229
    Ni = 1000
    R0 = 0.66
    B0 = 1.0
    R = np.linspace(0.05, 1.6, Ni)
    B = np.zeros(Ni)
    BIdeal = np.zeros(Ni)
    TF.Method = 'Filament'
    for i in range(Ni):
        Rvec = np.array([R[i], 0, 0])
        Bvec = TF.local(Rvec)
        B[i] = norm(Bvec)
    TF.Method = 'Simple'
    for i in range(Ni):
        Rvec = np.array([R[i], 0, 0])
        BIdeal[i] = norm(TF.local(Rvec))

    plt.plot(R, B)
    plt.plot(R, BIdeal, '-')
    plt.ylim(0, 1.2 * max(B))
    plt.title(r'Magnitude of B-field from TF Coils |B$_\phi(r)$|/|B$_\phi(R_o)$|')
    plt.xlabel('Radial Position [m]')
    plt.ylabel(r'Normalized B-Field Magnitude |B$_\phi$/B$_o$|')
    plt.legend(('Discrete Coils', r'Ideal B$\sim$1/R'))

# ------------------------------------------------------------------------------
# Plot Vertical Field in 1D at distance Z from the midplane


def PlotVF1D(VF):
    Ni = 1000
    R = np.linspace(0.15, 2.0, Ni)
    Z = [0.1, 0.2, 0.3]
    BNorm = np.zeros(Ni)
    BVert = np.zeros(Ni)
    BIdeal = np.zeros(Ni)
    Bvector = []
    Bmag = []
    for j in range(len(Z)):
        for i in range(Ni):
            Rvec = np.array([R[i], 0, Z[j]])
            Bvec = VF.local(Rvec)
            #BIdeal[i] = B0*R0/R[i]
            Bvector.append(Bvec)
            BNorm[i] = norm(Bvec) / VF.B0
            BVert[i] = Bvec[2] / VF.B0
#		B = B/B[1]
        plt.subplot(1, 2, 1)
        plt.plot(R, BNorm, label='Z = %0.1f' % Z[j])
        plt.subplot(1, 2, 2)
        plt.plot(R, BVert, label='Z = %0.1f' % Z[j])

#	plt.plot(R,BIdeal,'-');
    plt.subplot(1, 2, 1)
    plt.ylim(0, 1.2 * max(BNorm))
    plt.title(r'Magnitude of Vertical Field B$_z$($r$)')
    plt.xlabel('Radial Position [m]')
    plt.ylabel(r'Normalized B-Field Magnitude |B|/|B$_0$|')
    plt.legend(loc=3)

    plt.subplot(1, 2, 2)
    plt.ylim(1.2 * min(BVert), 1.2 * max(BVert))
    plt.title(r'Z-Component Vertical Field B$_z$($r$)')
    plt.xlabel('Radial Position [m]')
    plt.ylabel(r'Normalized Vertical Component of B-Field |B$_z$|/|B$_0$|')
    plt.legend(loc=3)

# ------------------------------------------------------------------------------
# Magnitude Plot of Toroidal Field in 2D


def PlotTF2D(TF):
    Nij = 100
    x = np.linspace(0, 1.4, Nij, float)
    y = np.linspace(0, 1.4, Nij, float)

    BMagMatrix = np.zeros((Nij, Nij), float)
    for i in range(Nij):
        print(i)
        for j in range(Nij):
            R = np.array([x[i], y[j]])
            B = TF.local(R)
            Mag = plt.norm(B)
            if (Mag < 2.5 and plt.norm(R) < 1.0) or (Mag < 0.75):
                BMagMatrix[i, j] = Mag
    plt.figure()
#	plt.pcolor(x,y,BMagMatrix); plt.colorbar
    plt.contour(x, y, BMagMatrix, 100)
    plt.title(r'Magnitude of B-field from TF Coils |B($r,\phi$)|/|B($R_o,\phi$)|')
    plt.xlabel('x [m]')
    plt.ylabel('y [m]')
    plt.colorbar()

# ------------------------------------------------------------------------------
# Magnitude Plot of Vertical Field in 2D


def PlotVF2D():
    Nij = 200
    r = np.linspace(0.2, 2.2, Nij)
    z = np.linspace(-1.0, 1.0, Nij)
    BMagMatrix = np.zeros((Nij, Nij), float)
    BZ0 = plt.norm(VF.local(np.array([0.1, 0, 0])))
    for i in range(Nij):
        print(i)
        for j in range(Nij):
            R = np.array([r[i], 0, z[j]])
            B = VF.local(R)
            Mag = plt.norm(B)
            if True:  # (Mag < 100.0):
                BMagMatrix[j, i] = Mag / BZ0  # log(Mag)
    plt.figure(4)
#	plt.pcolor(r,z,BMagMatrix);
    plt.contour(r, z, BMagMatrix, 120)
    plt.title(r'Magnitude of B-field from VF Coils: |B($r,z$)|/|B($0,0$)|')
    plt.xlabel('r [m]')
    plt.ylabel('z [m]')
    plt.colorbar()


# ------------------------------------------------------------------------------
# Vector plot of Vertical B-Field in 2D
def VectorVerticalField(VF):
    Nij = 25
#	r = np.linspace(0.2,2.2,Nij)
#	z = np.linspace(-1.0,1.0,Nij)
    r = np.linspace(1.0, 2.0, Nij)
    z = np.linspace(-0.5, 0.5, Nij)

    Br = []
    Bz = []
    R = []
    Z = []
    Mag = []
    for i in range(Nij):
        print(i)
        for j in range(Nij):
            R1 = np.array([r[i], 0, z[j]])
            B = VF.local(R1)
            R.append(r[i])
            Z.append(z[j])
            Br.append(B[0])
            Bz.append(B[2])
            Mag.append(np.sqrt(B[0]**2 + B[1]**2 + B[2]**2))
#	plt.figure()
    Q = plt.quiver(R, Z, Br, Bz, Mag, pivot='mid', scale=25,
                  width=0.005, cmap=mpl.cm.autumn)
    plt.title(r'Vertical Field Map (Poloidal View) |B$(r,z)$|/|B$(0,0)$|')
    plt.xlabel('r [m]')
    plt.ylabel('z [m]')
    plt.colorbar()
    plt.xlim(min(r), RInj[0])
    plt.ylim(min(z), max(z))
    return R, Z, Br, Bz, Mag

# ------------------------------------------------------------------------------
# Vector plot of Toroidal B-Field in 2D


def VectorToroidalField(TF):
    Nij = 25
#	x = np.linspace(0,1.4,Nij)
#	y = np.linspace(0,1.4,Nij)
#	x = np.linspace(1.0,1.4,Nij,float)
    x = np.linspace(0.9, RInj[0], Nij, float)
    y = np.linspace(-0.4, 0.4, Nij, float)

    X = []
    Y = []
    Bx = []
    By = []
    Mag = []
    for i in range(Nij):
        # print i
        for j in range(Nij):
            R = np.array([x[i], y[j], 0.0])
            B = TF.local(R)
            X.append(x[i])
            Y.append(y[j])
            Bx.append(B[0])
            By.append(B[1])
            MAG = np.sqrt(B[0]**2 + B[1]**2 + B[2]**2)
            if (MAG < 0.8 and plt.norm(R) < 1.0) or (MAG < 0.8):
                Mag.append(MAG)
            else:
                Mag.append(nan)
    Q = plt.quiver(X, Y, Bx, By, np.array(Mag), pivot='mid', scale=10,
                  width=0.005, cmap=mpl.cm.winter)  # ,cmap=mpl.cm.winter
#	plt.title(r'Toroidal B($r,\phi$)-field (uniform in $z$)')
    plt.title(r'Toroidal Field Map (Top View) B($r,\phi$)/|B$(R_o,\phi)$|')
    plt.xlabel('x [m]')
    plt.ylabel('y [m]')
    cb = plt.colorbar()
#	plt.xlim(min(x),max(x));plt.ylim(min(y),max(y))
    plt.xlim(min(x), RInj[0])
    plt.ylim(min(y), max(y))
    return X, Y, Bx, By, Mag


#RInj = [1.798, -0.052, 0.243]

# ------------------------------------------------------------------------------
# Draw Beam
alpha = 12.6 / 180.0 * pi
beta = 8.0 / 180.0 * pi
#Rinjection = [1.798, 0.052, 0.243]
xTop = [RInj[0], 0.6]
yTop = [RInj[1], RInj[1] + (xTop[0] - xTop[1]) * np.tan(beta)]

rTop = [2.0, 1.798, 0.6]
zTop = [0.243 + (rTop[0] - rTop[1]) * np.tan(alpha), 0.243,
        0.243 - (rTop[1] - rTop[2]) * np.tan(alpha)]


if True:
    plt.figure(1)
    X, Y, Bx, By, Mag = VectorToroidalField(TF)
    plt.plot(xTop, yTop, 'r', linewidth=2)

if True:
    plt.figure(2, figsize=(16, 6))

    plt.subplot(1, 2, 1)
    X, Y, Bx, By, Mag = VectorToroidalField(TF)
    plt.subplot(1, 2, 1)
    plt.plot(xTop, yTop, 'r', linewidth=2)
    plt.xlabel('x [m]')
    plt.ylabel('y [m]')
    plt.title('Toroidal Field')

    plt.subplot(1, 2, 2)
    R, Z, Br, Bz, Mag = VectorVerticalField(VF)
    plt.subplot(1, 2, 2)
    plt.plot(rTop, zTop, 'r', linewidth=2)
    plt.xlabel('R [m]')
    plt.ylabel('Z [m]')
    plt.title('Vertical Field')

if True:
    plt.figure(3)
    PlotTF1D(TF)

if True:
    plt.figure(4)
    PlotVF1D(VF)

plt.show()
