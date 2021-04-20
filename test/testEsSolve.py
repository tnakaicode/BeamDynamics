import numpy as np
import math
import os
import sys
sys.path.append(os.path.join('../'))

from solvers.esSolve import electricFieldAtPoint, potentialToElectricField
from solvers.esSolve import laplace, laplace1D, laplace2D, laplace3D
from solvers.neumann import neumann as neuBC
from solvers.dirichlet import dirichlet as dirBC

# 1D grid
NX_1D = 10  # gauss seidel is the one that's breaking? with jacobi only, it will break when it goes to 50
LX_1D = 1.2
DX_1D = LX_1D / NX_1D
X0_1D = 1.0

# 2D grid
NX_2D = 6
NY_2D = 7
LX_2D = 0.01 * 1.25
LY_2D = 0.01 * 2.3
DX_2D = LX_2D / NX_2D
DY_2D = LY_2D / NY_2D
X0_2D = 1.0
Y0_2D = 2.0

# 3D grid
NX_3D = 5
NY_3D = 6
NZ_3D = 7
LX_3D = 0.01 * 1.25
LY_3D = 0.01 * 2.3
LZ_3D = 0.01 * 1.87
DX_3D = LX_3D / NX_3D
DY_3D = LY_3D / NY_3D
DZ_3D = LZ_3D / NZ_3D
X0_3D = 1.0
Y0_3D = 2.0
Z0_3D = 3.0

# 1D test cases


def V1_1D(x):
    return 1.0 * x + 2.0


def E1_1D(x):
    return -np.array([1.0])

# 2D test cases


def V1_2D(x, y):
    return 0.5 * (pow(x, 2.0) - pow(y, 2.0))


def E1_2D(x, y):
    return -np.array([x, -y])


def V2_2D(x, y):
    return 2.0 * x + 0.5 * y


def E2_2D(x, y):
    return -np.array([2.0, 0.5])


def V3_2D(x, y):
    return 2.0


def E3_2D(x, y):
    return -np.array([0.0, 0.0])


def V4_2D(x, y):
    a = 0.1 * 2.0 * math.pi / LX_2D
    return (math.cos(a * x) + math.sin(a * x)) * (math.cosh(a * y) + math.sinh(a * y))


def E4_2D(x, y):
    a = 0.1 * 2.0 * math.pi / LX_2D
    return -a * np.array([(math.cos(a * x) - math.sin(a * x)) * (math.cosh(a * y) + math.sinh(a * y)),
                          (math.cos(a * x) + math.sin(a * x)) * (math.cosh(a * y) + math.sinh(a * y))])


def V5_2D(x, y):
    return V1_2D(x, y) + V2_2D(x, y) + V3_2D(x, y)


def E5_2D(x, y):
    return np.array([E1_2D(x, y)[0] + E2_2D(x, y)[0] + E3_2D(x, y)[0],
                     E1_2D(x, y)[1] + E2_2D(x, y)[1] + E3_2D(x, y)[1]])

# 3D test cases


def V1_3D(x, y, z):
    return 0.5 * V5_2D(x, y) + 2.0 * z


def E1_3D(x, y, z):
    return np.array([0.5 * E5_2D(x, y)[0], 0.5 * E5_2D(x, y)[1], -2.0])


def V2_3D(x, y, z):
    return 0.5 * V5_2D(x, z) + 2.0 * y


def E2_3D(x, y, z):
    return np.array([0.5 * E5_2D(x, z)[0], -2.0, 0.5 * E5_2D(x, z)[1]])


def V3_3D(x, y, z):
    return 0.5 * V5_2D(y, z) + 2.0 * x


def E3_3D(x, y, z):
    return np.array([-2.0, 0.5 * E5_2D(y, z)[0], 0.5 * E5_2D(y, z)[1]])


def V4_3D(x, y, z):
    a = 0.1 * 2.0 * math.pi / max(LX_3D, LY_3D)
# should these be using np versions?
    return math.exp(a * x) * math.exp(a * y) * math.sin(math.sqrt(2.0) * a * z)


def E4_3D(x, y, z):
    a = 0.1 * 2.0 * math.pi / max(LX_3D, LY_3D)
# should these be using np versions?
    return -np.array([a * math.exp(a * x) * math.exp(a * y) * math.sin(math.sqrt(2.0) * a * z),
                      a * math.exp(a * x) * math.exp(a * y) *
                      math.sin(math.sqrt(2.0) * a * z),
                      math.sqrt(2.0) * a * math.exp(a * x) * math.exp(a * y) * math.cos(math.sqrt(2.0) * a * z)])


def testElectricFieldAtPoint():
    absTol = 1.0e-5
    relTol = 0.0

    def test(array1, array2):
        assert np.allclose(array1, array2, relTol, absTol)

    def test1D():
        def E(x):
            return math.exp(x)

        E_grid = np.zeros((NX_1D + 1, 1))

        for i in range(NX_1D + 1):
            E_grid[i] = E(X0_1D + DX_1D * i)

        numPointsToTestBetweenGridPoints = 10

        for i in range(numPointsToTestBetweenGridPoints * NX_1D + 1):
            Xi = X0_1D + DX_1D * i / numPointsToTestBetweenGridPoints

        test(electricFieldAtPoint(
            E_grid, [DX_1D], [X0_1D], [Xi]), E(Xi))

    def test2D():
        def E(x, y):
            return math.exp(x + y)

        E_grid = np.zeros((NX_2D + 1, NY_2D + 1, 2))

        for i in range(NX_2D + 1):
            for j in range(NY_2D + 1):
                E_grid[i][j] = E(X0_2D + DX_2D * i, Y0_2D + DY_2D * j)

        numPointsToTestBetweenGridPoints_x = 10
        numPointsToTestBetweenGridPoints_y = 10

        for i in range(numPointsToTestBetweenGridPoints_x * NX_2D + 1):
            for j in range(numPointsToTestBetweenGridPoints_y * NY_2D + 1):
                Xi = X0_2D + DX_2D * i / numPointsToTestBetweenGridPoints_x
                Yi = Y0_2D + DY_2D * j / numPointsToTestBetweenGridPoints_y

        test(electricFieldAtPoint(
            E_grid, [DX_2D, DY_2D], [X0_2D, Y0_2D], [Xi, Yi]), E(Xi, Yi))

    def test3D():
        def E(x, y, z):
            return math.exp(x + y + z)

        E_grid = np.zeros((NX_3D + 1, NY_3D + 1, NZ_3D, 3))

        for i in range(NX_3D + 1):
            for j in range(NY_3D + 1):
                for k in range(NY_3D + 1):
                    E_grid[i][j][k] = E(X0_3D + DX_3D * i,
                                        Y0_3D + DY_3D * j, Z0_3D + DZ_3D * k)

        numPointsToTestBetweenGridPoints_x = 10
        numPointsToTestBetweenGridPoints_y = 10
        numPointsToTestBetweenGridPoints_z = 10

        for i in range(numPointsToTestBetweenGridPoints_x * NX_3D + 1):
            for j in range(numPointsToTestBetweenGridPoints_y * NY_3D + 1):
                for k in range(numPointsToTestBetweenGridPoints_z * NY_3D + 1):
                    Xi = X0_3D + DX_3D * i / numPointsToTestBetweenGridPoints_x
                    Yi = Y0_3D + DY_3D * j / numPointsToTestBetweenGridPoints_y
                    Zi = Z0_3D + DZ_3D * k / numPointsToTestBetweenGridPoints_z

        test(electricFieldAtPoint(E_grid, [DX_3D, DY_3D, DZ_3D], [
             X0_3D, Y0_3D, Z0_3D], [Xi, Yi, Zi]), E(Xi, Yi, Zi))

    test1D()
    test2D()
    test3D()


def testPotentialToElectricField():
    absTol = 1.0e-2
    relTol = 0.0

    def test(array1, array2):
        assert np.allclose(array1, array2, relTol, absTol)

    def test1D(potential, field):
        potAccept1D = np.zeros(NX_1D + 1)
        fieldAccept1D = np.zeros((NX_1D + 1, 1))
        for i in range(NX_1D + 1):
            potAccept1D[i] = potential(X0_1D + DX_1D * i)
            fieldAccept1D[i] = field(X0_1D + DX_1D * i)

        test(potentialToElectricField(potAccept1D, [DX_1D]),
             fieldAccept1D)

    def test2D(potential, field):
        potAccept2D = np.zeros((NX_2D + 1, NY_2D + 1))
        fieldAccept2D = np.zeros((NX_2D + 1, NY_2D + 1, 2))
        for i, j in np.ndindex(potAccept2D.shape):
            potAccept2D[i, j] = potential(X0_2D + DX_2D * i, Y0_2D + DY_2D * j)
            fieldAccept2D[i, j] = field(X0_2D + DX_2D * i, Y0_2D + DY_2D * j)

        test(potentialToElectricField(potAccept2D, [DX_2D, DY_2D]),
             fieldAccept2D)

    def test3D(potential, field):
        potAccept3D = np.zeros((NX_3D + 1, NY_3D + 1, NZ_3D + 1))
        fieldAccept3D = np.zeros((NX_3D + 1, NY_3D + 1, NZ_3D + 1, 3))
        for i, j, k in np.ndindex(potAccept3D.shape):
            # consider using np.fromfunction here
            potAccept3D[i, j, k] = potential(
                X0_3D + DX_3D * i, Y0_3D + DY_3D * j, Z0_3D + DZ_3D * k)
            fieldAccept3D[i, j, k] = field(
                X0_3D + DX_3D * i, Y0_3D + DY_3D * j, Z0_3D + DZ_3D * k)

        test(potentialToElectricField(potAccept3D, [DX_3D, DY_3D, DZ_3D]),
             fieldAccept3D)

    # Run the tests
    test1D(V1_1D, E1_1D)

    # , [V4_2D,E4_2D]]:
    for testFuncs in [[V1_2D, E1_2D], [V2_2D, E2_2D], [V3_2D, E3_2D], [V5_2D, E5_2D]]:
        test2D(testFuncs[0], testFuncs[1])

    # , [V4_3D,E4_3D]]:
    for testFuncs in [[V1_3D, E1_3D], [V2_3D, E2_3D], [V3_3D, E3_3D]]:
        test3D(testFuncs[0], testFuncs[1])


def test_laplace():
    absTol = 1.0e-3
    relTol = 0.0

    def test(array1, array2):
        assert np.allclose(array1, array2, relTol, absTol)

    def test1D(potential, field):
        potAccept1D = np.zeros(NX_1D + 1)
        fieldAccept1D = np.zeros((NX_1D + 1, 1))
        for i in range(NX_1D + 1):
            potAccept1D[i] = potential(X0_1D + DX_1D * i)
            fieldAccept1D[i] = field(X0_1D + DX_1D * i)

        # Boundary conditions
        V0_1D = dirBC(potAccept1D[0])
        VN_1D = dirBC(potAccept1D[NX_1D])

        E0_1D = neuBC(fieldAccept1D[0])
        EN_1D = neuBC(fieldAccept1D[NX_1D])

        for BC0_1D in [V0_1D, E0_1D]:
            for BCN_1D in [VN_1D]:  # ,EN_1D]:
                for testType in ["direct", "jacobi", "gaussSeidel"]:
                    for cythonType in [True, False]:
                        if testType == "direct":
                            potCalculated1D = laplace1D(NX_1D, DX_1D, BC0_1D, BCN_1D,
                                                        testType, cythonType)
                        else:
                            potCalculated1D = laplace1D(NX_1D, DX_1D, BC0_1D, BCN_1D,
                                                        testType, relTol, absTol, cythonType)

                        test(potCalculated1D, potAccept1D)

    def test2D(potential, field):
        potAccept2D = np.zeros((NX_2D + 1, NY_2D + 1))
        fieldAccept2D = np.zeros((NX_2D + 1, NY_2D + 1, 2))
        for i, j in np.ndindex(potAccept2D.shape):
            potAccept2D[i, j] = potential(X0_2D + DX_2D * i, Y0_2D + DY_2D * j)
            fieldAccept2D[i, j] = field(X0_2D + DX_2D * i, Y0_2D + DY_2D * j)

        # Boundary conditions
        V0x_2D = dirBC(potAccept2D[0, :])
        VNx_2D = dirBC(potAccept2D[NX_2D, :])
        V0y_2D = dirBC(potAccept2D[:, 0])
        VNy_2D = dirBC(potAccept2D[:, NY_2D])

        E0x_2D = neuBC(fieldAccept2D[0, :, 0])
        ENx_2D = neuBC(fieldAccept2D[NX_2D, :, 0])
        E0y_2D = neuBC(fieldAccept2D[:, 0, 1])
        ENy_2D = neuBC(fieldAccept2D[:, NY_2D, 1])

        for BC0x_2D in [V0x_2D, E0x_2D]:
            for BCNx_2D in [VNx_2D, ENx_2D]:
                for BC0y_2D in [V0y_2D, E0y_2D]:
                    for BCNy_2D in [VNy_2D]:  # ,ENy_2D]:
                        for testType in ["direct", "jacobi", "gaussSeidel"]:
                            for cythonType in [True, False]:
                                if testType == "direct":
                                    potCalculated2D = laplace2D(NX_2D, DX_2D, BC0x_2D, BCNx_2D,
                                                                NY_2D, DY_2D, BC0y_2D, BCNy_2D,
                                                                testType, cythonType)
                                else:
                                    potCalculated2D = laplace2D(NX_2D, DX_2D, BC0x_2D, BCNx_2D,
                                                                NY_2D, DY_2D, BC0y_2D, BCNy_2D,
                                                                testType, relTol, absTol, cythonType)

                                test(potCalculated2D, potAccept2D)

    def test3D(potential, field):
        potAccept3D = np.zeros((NX_3D + 1, NY_3D + 1, NZ_3D + 1))
        fieldAccept3D = np.zeros((NX_3D + 1, NY_3D + 1, NZ_3D + 1, 3))
        for i, j, k in np.ndindex(potAccept3D.shape):
            # consider using np.fromfunction here
            potAccept3D[i, j, k] = potential(
                X0_3D + DX_3D * i, Y0_3D + DY_3D * j, Z0_3D + DZ_3D * k)
            fieldAccept3D[i, j, k] = field(
                X0_3D + DX_3D * i, Y0_3D + DY_3D * j, Z0_3D + DZ_3D * k)

        # Boundary conditions
        V0x_3D = dirBC(potAccept3D[0, :, :])
        VNx_3D = dirBC(potAccept3D[NX_3D, :, :])
        V0y_3D = dirBC(potAccept3D[:, 0, :])
        VNy_3D = dirBC(potAccept3D[:, NY_3D, :])
        V0z_3D = dirBC(potAccept3D[:, :, 0])
        VNz_3D = dirBC(potAccept3D[:, :, NZ_3D])

        E0x_3D = neuBC(fieldAccept3D[0, :, :, 0])
        ENx_3D = neuBC(fieldAccept3D[NX_3D, :, :, 0])
        E0y_3D = neuBC(fieldAccept3D[:, 0, :, 1])
        ENy_3D = neuBC(fieldAccept3D[:, NY_3D, :, 1])
        E0z_3D = neuBC(fieldAccept3D[:, :, 0, 2])
        ENz_3D = neuBC(fieldAccept3D[:, :, NZ_3D, 2])

        for BC0x_3D in [V0x_3D, E0x_3D]:
            for BCNx_3D in [VNx_3D, ENx_3D]:
                for BC0y_3D in [V0y_3D, E0y_3D]:
                    for BCNy_3D in [VNy_3D, ENy_3D]:
                        for BC0z_3D in [V0z_3D, E0z_3D]:
                            for BCNz_3D in [VNz_3D]:  # ,ENz_3D]:
                                for testType in ["direct", "jacobi", "gaussSeidel"]:
                                    for cythonType in [True, False]:
                                        if testType == "direct":
                                            potCalculated3D = laplace3D(NX_3D, DX_3D, BC0x_3D, BCNx_3D,
                                                                        NY_3D, DY_3D, BC0y_3D, BCNy_3D,
                                                                        NZ_3D, DZ_3D, BC0z_3D, BCNz_3D,
                                                                        testType, cythonType)
                                        else:
                                            potCalculated3D = laplace3D(NX_3D, DX_3D, BC0x_3D, BCNx_3D,
                                                                        NY_3D, DY_3D, BC0y_3D, BCNy_3D,
                                                                        NZ_3D, DZ_3D, BC0z_3D, BCNz_3D,
                                                                        testType, relTol, absTol, cythonType)

                                        test(potCalculated3D, potAccept3D)

    # Run the tests
    test1D(V1_1D, E1_1D)

    # , [V4_2D,E4_2D]]:
    for testFuncs in [[V1_2D, E1_2D], [V2_2D, E2_2D], [V3_2D, E3_2D], [V5_2D, E5_2D]]:
        test2D(testFuncs[0], testFuncs[1])

    # , [V4_3D,E4_3D]]:
    for testFuncs in [[V1_3D, E1_3D], [V2_3D, E2_3D], [V3_3D, E3_3D]]:
        test3D(testFuncs[0], testFuncs[1])
