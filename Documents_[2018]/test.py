import numpy as np
import fmodpy

test = fmodpy.fimport("test.f90")

np.random.seed(0)
a = np.array(list(range(10)), dtype=np.float64)
# i = np.array([1,2,4], dtype=np.int64)
a[[3,5,7]] = 0
print(a)
test.nonzero(a)

print()
a = np.random.random(shape=(3,))
b = np.random.random(shape=(3,4))
c = test.vec_mat(a, b)
print(a)
print(b)
print(c)
