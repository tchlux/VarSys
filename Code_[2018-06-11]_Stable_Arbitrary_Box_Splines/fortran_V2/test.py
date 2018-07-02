import numpy as np

# A matrix of direction (column) vectors
mat = np.array([[1.,0., 1., 1.],
                [0.,1.,-1., 1.]])
print()
print(mat)
print()

# Numerically unstable solution (but faster to compute)
min_norm = np.linalg.solve(np.matmul(mat,mat.T), mat)
print(min_norm)
print()

# Numerically stable solution (but slower to compute)
min_norm, residuals, rank, singular_vals = np.linalg.lstsq(
    mat.T, np.identity(mat.shape[1]), rcond=None)
print(min_norm)
print()

print('-'*70)

mat = min_norm[:,1:]
print()
print(mat)
print()

# Numerically unstable solution (but faster to compute)
min_norm = np.linalg.solve(np.matmul(mat,mat.T), mat)
print(min_norm)
print()

# Numerically stable solution (but slower to compute)
min_norm, residuals, rank, singular_vals = np.linalg.lstsq(
    mat.T, np.identity(mat.shape[1]), rcond=None)
print(min_norm)
print()

import fmodpy
fmodpy.wrap("min_norm.f90", module_link_args=["-lblas","-llapack","-lgfortran"])
import min_norm
min_norm.matrix_minimum_norm(mat)
