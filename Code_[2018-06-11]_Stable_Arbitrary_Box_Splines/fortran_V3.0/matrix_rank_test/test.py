import fmodpy

code = fmodpy.fimport("code.f90", verbose=True, 
                      module_link_args=["-lgfortran", "-lblas", "-llapack"])

# Generate some test data
import numpy as np
# Make the arrays f-contiguous
a = np.asarray([
    [1, 0, 0, 0],
    [0, 1, 0, 0],
    [1, 1, 1, 1],
    [1, 1, 1, 1]
], order='F', dtype=np.float64)

print(a)
print(code.matrix_rank(a))
# s = np.linalg.svd(a, compute_uv=False)
# print(s)
# tol = s.max() * max(a.shape) * np.eps
# print(np.linalg.matrix_rank(a))
