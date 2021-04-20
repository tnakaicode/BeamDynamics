import cython
import numpy as np
cimport numpy as np

# Solves A x = B iteratively
@cython.boundscheck(False)
@cython.wraparound(False)
def iterate(np.ndarray[double, ndim=2] A, np.ndarray[double, ndim=1] B, int solType, double relTol, double absTol):
  cdef unsigned int maxIterations = 1000000
  cdef np.ndarray[double, ndim=1] x = np.copy(B)
  cdef np.ndarray[double, ndim=1] oldX
  cdef Py_ssize_t i
  cdef Py_ssize_t iteration
  for iteration in range(maxIterations):
    oldX = np.copy(x)
    for i in range(x.shape[0]):
      # should really first try to rearrange A,B to make A[i,i] != 0 if possible
      if A[i,i] != 0.0:
        if (solType == 1):
          x[i] = (1.0/A[i][i])*(B[i] - np.dot(A[i,:i],oldX[:i]) - np.dot(A[i,i+1:], oldX[i+1:]))
        elif (solType == 2):
          x[i] = (1.0/A[i][i])*(B[i] - np.dot(A[i,:i],x[:i]) - np.dot(A[i,i+1:], oldX[i+1:]))
    if np.allclose(oldX,x,relTol,absTol):
      break
  return x

# Solves A x = B using Jacobi iteration
def jacobi(np.ndarray[double, ndim=2] A, np.ndarray[double, ndim=1] B, double relTol, double absTol):
  return iterate(A,B,1,relTol/1000.0,absTol/1000.0)

# Solves A x = B using Gauss Seidel iteration
def gaussSeidel(np.ndarray[double, ndim=2] A, np.ndarray[double, ndim=1] B, double relTol, double absTol):
  return iterate(A,B,2,relTol/100.0,absTol/100.0)

# Solves A x = B directly
def direct(np.ndarray[double, ndim=2] A, np.ndarray[double, ndim=1] B):
  return np.linalg.solve(A,B)
