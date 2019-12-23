import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sp
import sys
import os
sys.path.append(os.path.join('../'))
from lib.BeamDynamicsTools.Boundary import Boundary
from lib.BeamDynamicsTools.Bfield import Bfield, BfieldTF, BfieldVF
from lib.BeamDynamicsTools.Trajectory import Trajectory
from lib.BeamDynamicsTools.Beam import Beam
from lib.BeamDynamicsTools.Ellipse import Ellipse
import numpy as np
from matplotlib.pyplot import *

Path0 = '../output/sigma_trace3d/igmaFinal/'
FileList = ['SigmaFinal_I_0.txt', 'SigmaFinal_I_1600.txt',
            'SigmaFinal_I_3120.txt', 'SigmaFinal_I_4450.txt']
FileListT3D = ['Trace3DSigma_I_0.txt', 'Trace3DSigma_I_1600NG.txt',
               'Trace3DSigma_I_3120NG.txt', 'Trace3DSigma_I_4450NG.txt']
Color = ['b', 'g', 'r', 'c']
LIM = 12

# ======== Import Sigma Matrices ==============================================
Sigma = []
Ellipse = []
SigmaX = []
SigmaY = []
SigmaT3D = []
EllipseT3D = []
SigmaXT3D = []
SigmaYT3D = []

SigmaLC = np.loadtxt(Path0 + 'Trace3DSigma_I_0LC.txt')
SigmaInj = np.matrix(np.loadtxt(Path0 + 'SigmaInj.txt'))
ellipseInj = Ellipse(SigmaInj)

for i in range(len(FileList)):
    Sigma.append(np.matrix(np.loadtxt(Path0 + FileList[i])))
    sp.ellipse.append(Ellipse(Sigma[-1]))

for i in range(len(FileListT3D)):
    SigmaT3D.append(np.matrix(np.loadtxt(Path0 + FileListT3D[i])))
    EllipseT3D.append(Ellipse(SigmaT3D[-1]))

# ======== Calculate Mismatch =================================================

M0000 = sp.ellipse[0].MismatchFactor(EllipseT3D[0])
M1600 = sp.ellipse[1].MismatchFactor(EllipseT3D[1])
M3120 = sp.ellipse[2].MismatchFactor(EllipseT3D[2])
M4500 = sp.ellipse[3].MismatchFactor(EllipseT3D[3])

M = []
for i in [0, 1, 2, 3]:
    M.append(EllipseT3D[i].MismatchFactor(sp.ellipse[i]))

# ======== All in one plot ====================================================
if False:
    plt.figure(figsize=(12, 4))
    C1 = 1.0  # /np.sqrt(5.0)#0.9
    for i in [0, 1]:
        subplot(1, 3, 1, aspect='equal')
        sp.ellipse[i].PlotXX1(Mod=Color[i])
        EllipseT3D[i].PlotXX1(Mod='--' + Color[i])
        subplot(1, 3, 2, aspect='equal')
        sp.ellipse[i].PlotYY1(Mod=Color[i])
        EllipseT3D[i].PlotYY1(Mod='--' + Color[i])
        subplot(1, 3, 3, aspect='equal')
        sp.ellipse[i].PlotXY(Mod=Color[i])
        EllipseT3D[i].PlotXY(Mod='--' + Color[i], Rotate=False)
        xlim(-LIM, LIM)
        ylim(-LIM, LIM)
        title(FileList[i])

# ======== All on separate plots ============================================
Title = 'Dynamics Modeling in C-Mod'
BList = ['0.0000 T', '0.0582 T', '0.1135 T', '0.1618 T']

if True:
    #	for i in [0,1,2,3]:
    for i in [2]:
        plt.figure(figsize=(10, 10))
#		title(FileList[i])
        subplot(2, 2, 1, aspect='equal')
        ellipseInj.PlotXX1(Mod=':r')
        sp.ellipse[i].PlotXX1(Mod=Color[0])
        EllipseT3D[i].PlotXX1(Mod='--' + Color[1])
        text(0, 0, 'M=%0.4f' % M[i][0], va='center',
             ha='center', color='r', size=16)
        title('Transverse Phase Plane (horizontal)')
        legend(('Initial Beam', 'AIMS Code', 'TRACE3D'), loc=2)

        subplot(2, 2, 2, aspect='equal')
        ellipseInj.PlotYY1(Mod=':r')
        sp.ellipse[i].PlotYY1(Mod=Color[0])
        EllipseT3D[i].PlotYY1(Mod='--' + Color[1])
        text(0, 0, 'M=%0.4f' % M[i][1], va='center',
             ha='center', color='r', size=16)
        title('Transverse Phase Plane (Vertical)')

        subplot(2, 2, 3, aspect='equal')
        ellipseInj.PlotZZ1(Mod=':r')
        sp.ellipse[i].PlotZZ1(Mod=Color[0])
        EllipseT3D[i].PlotZZ1(Mod='--' + Color[1])
        text(0, 0, 'M=%0.4f' % M[i][2], va='center',
             ha='center', color='r', size=16)
        title('Longitudinal Phase Plane')

        subplot(2, 2, 4, aspect='equal')
        ellipseInj.PlotXY(Mod=':r')
        sp.ellipse[i].PlotXY(Mod=Color[0], Rotate=False)
        EllipseT3D[i].PlotXY(Mod='--' + Color[1])
        text(0, 0, 'M=%0.4f' % M[i][3], va='center',
             ha='center', color='r', size=16)
        title('Transverse Projection')

        E = sp.ellipse[i]
        np.savetxt('Maxima.txt', [E.WidthX, E.WidthY, E.WidthZ, E.DivergenceX,
                               E.DivergenceY, E.DivergenceZ, M[i][0], M[i][1], M[i][2], M[i][3]])
        suptitle(Title + r' (B$_\phi$= ' + BList[i] + ')', size=16)

# ======== 90 degree bend test ============================================

if True:
    Sigma90 = np.loadtxt(Path0 + 'SigmaBend90.txt')
    Sigma90T3D = np.loadtxt(Path0 + 'Trace3DSigmaBend90.txt')

    E90 = Ellipse(Sigma90)
    E90T3D = Ellipse(Sigma90T3D)
    M90 = E90.MismatchFactor(E90T3D)

    plt.figure(figsize=(10, 10))
    FormatT = dict(va='center', ha='center', color='r', size=16)

    subplot(2, 2, 1, aspect='equal')
    ellipseInj.PlotXX1(Mod=':r')
    p1 = E90.PlotXX1(Mod=Color[0])
    p2 = E90T3D.PlotXX1(Mod='--' + Color[1])
    text(0, 0, 'M=%0.4f' % M90[0], FormatT)
    title('Transverse Phase Plane (horizontal)')
    legend(('Initial Beam', 'AIMS Code', 'TRACE3D'), loc=2)

    subplot(2, 2, 2, aspect='equal')
    ellipseInj.PlotYY1(Mod=':r')
    E90.PlotYY1(Mod=Color[0])
    E90T3D.PlotYY1(Mod='--' + Color[1])
    text(0, 0, 'M=%0.4f' % M90[1], FormatT)
    title('Transverse Phase Plane (Vertical)')

    subplot(2, 2, 3, aspect='equal')
    ellipseInj.PlotZZ1(Mod=':r')
    E90.PlotZZ1(Mod=Color[0])
    E90T3D.PlotZZ1(Mod='--' + Color[1])
    text(0, 0, 'M=%0.4f' % M90[2], FormatT)
    title('Longitudinal Phase Plane')

    subplot(2, 2, 4, aspect='equal')
    ellipseInj.PlotXY(Mod=':r')
    E90.PlotXY(Mod=Color[0], Rotate=False)
    E90T3D.PlotXY(Mod='--' + Color[1])
    text(0, 0, 'M=%0.4f' % M90[3], FormatT)
    title('Transverse Projection')

    suptitle(r'Dynamics Comparison (90$^o$Bend, Uniform B, R$_c$ = 1 m)', size=16)

    print(E90.SpatialWidth())
    print(E90.AngularWidth())
    print(M90)

# ======== Space Charge Effects ============================================

if False:
    plt.figure(figsize=(8, 8))
    E0 = Ellipse(Sigma[0])
    E1 = Ellipse(SigmaLC)
    M = E0.MismatchFactor(E1, Type=1)

    subplot(2, 2, 1, aspect='equal')
    E0.PlotXX1()
    E1.PlotXX1()
    text(0, 0, 'M=%0.4f' % M[1], va='center', ha='center', color='r', size=16)
    legend((r'1.000 mA', '0.001 mA'), loc=2)

    subplot(2, 2, 2, aspect='equal')
    E0.PlotYY1()
    E1.PlotYY1()
    text(0, 0, 'M=%0.4f' % M[1], va='center', ha='center', color='r', size=16)

    subplot(2, 2, 3, aspect='equal')
    E0.PlotZZ1()
    E1.PlotZZ1()
    text(0, 0, 'M=%0.4f' % M[2], va='center', ha='center', color='r', size=16)

    subplot(2, 2, 4, aspect='equal')
    E0.PlotXY()
    E1.PlotXY()
    text(0, 0, 'M=%0.4f' % M[3], va='center', ha='center', color='r', size=16)

    suptitle(r'Space Charge Effects: 1mA versus 1$\mu$A', size=16)

# ======== Plot All Projections ============================================

if True:
    plt.figure()
    for i in [0, 1, 2, 3]:
        sp.ellipse[i].PlotXY()

show()
