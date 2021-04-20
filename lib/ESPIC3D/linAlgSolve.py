import numpy as np

# Solves A x = B iteratively


def iterate(A, B, solType, relTol, absTol):
    maxIterations = 1000000
    x = np.copy(B)
    for iteration in range(maxIterations):
        oldX = np.copy(x)
        for i in range(x.shape[0]):
            # should really first try to rearrange A,B to make A[i,i] != 0 if possible
            if A[i, i] != 0.0:
                if (solType == "jacobi"):
                    x[i] = (1.0 / A[i][i]) * (B[i] - np.dot(A[i, :i],
                                                            oldX[:i]) - np.dot(A[i, i + 1:], oldX[i + 1:]))
                elif (solType == "gaussSeidel"):
                    x[i] = (1.0 / A[i][i]) * (B[i] - np.dot(A[i, :i],
                                                            x[:i]) - np.dot(A[i, i + 1:], oldX[i + 1:]))
        if np.allclose(oldX, x, relTol, absTol):
            break
    return x

# Solves A x = B using Jacobi iteration


def jacobi(A, B, relTol, absTol):
    return iterate(A, B, "jacobi", relTol / 1000.0, absTol / 1000.0)

# Solves A x = B using Gauss Seidel iteration


def gaussSeidel(A, B, relTol, absTol):
    return iterate(A, B, "gaussSeidel", relTol / 100.0, absTol / 100.0)

# Solves A x = B directly


def direct(A, B):
    return np.linalg.solve(A, B)
