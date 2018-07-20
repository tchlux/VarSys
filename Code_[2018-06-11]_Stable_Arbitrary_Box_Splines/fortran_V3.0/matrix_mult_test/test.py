import fmodpy

code = fmodpy.fimport("code.f90", verbose=True, 
                      module_link_args=["-lgfortran", "-lblas", "-llapack"])

# Generate some test data
import numpy as np
a = np.linspace(0,7,8).reshape((2,4))
b = np.linspace(0,1,2).reshape((2,1))
c = np.ones((4,1))
# Make the arrays f-contiguous
a = np.asarray(a, order='F')
b = np.asarray(b, order='F')
c = np.asarray(c, order='F')
# Printout initial samples
print()
print(a)
print(b)
print(c)
print()
# Execute code
code.mm(a,b,c)
# Print output
print(c)
