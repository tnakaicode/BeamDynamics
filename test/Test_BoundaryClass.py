import sys
sys.path.append('../lib/BeamDynamicsTools/')
from Boundary import *
import numpy as np
import pylab as pl
import timeit

# Import poloidal boundary points
Rb = np.loadtxt('../data/CmodCoordinatesRZ.dat', usecols=[0])
Zb = np.loadtxt('../data/CmodCoordinatesRZ.dat', usecols=[1])

# ------------------------------------------------------------------------------
# Generate vessel boundary
Vessel = Boundary(Rb, Zb)

# ------------------------------------------------------------------------------
plt.figure(1)
Vessel.Border()

plt.figure(2)
Vessel.Plot2D(2)

# ------------------------------------------------------------------------------
# 3D plot of vessel boundary
ax = Vessel.Figure3D(3)
Vessel.Plot3D(ax)


# ===============================================================================
# Test Case for in-out detection algorithm
# ===============================================================================

plt.figure(4)
Xrand = np.array([0.0, 0.0])
Ni = 10000
Xrange = [0.2, 1.2]
Zrange = [0.8, -0.8]

# ------------------------------------------------------------------------------
# Generate Ni Random X,Z points
Xrand = random.rand(Ni) * (Xrange[1] - Xrange[0]) + Xrange[0]
Zrand = random.rand(Ni) * (Zrange[1] - Zrange[0]) + Zrange[0]

# ------------------------------------------------------------------------------
# Test Each Point to determine if is in the boundary

IN = []
start = timeit.default_timer()
for i in range(Ni):
    r = [Xrand[i], 0.0, Zrand[i]]
    IN.append(Vessel.InBoundary(r))
stop = timeit.default_timer()
TIME = (stop - start) / Ni * 1000

InX = []
InY = []
OutX = []
OutZ = []
for i in range(Ni):
    if IN[i]:
        InX.append(Xrand[i])
        InY.append(Zrand[i])
    else:
        OutX.append(Xrand[i])
        OutZ.append(Zrand[i])

# ------------------------------------------------------------------------------
# Plot Results

plt.plot(OutX, OutZ, '.r')
plt.plot(InX, InY, '.g')
plt.legend(('Out', 'In'))
plt.title(
    'In-Out Boundary Detection: Time %0.4f ms/test (N$_{test}$ = %0.0f)' % (TIME, Ni))
Vessel.Border()


plt.show()
