

# Import the splines Fortran code.
import fmodpy
# Import the TOMS code.
bspline = fmodpy.fimport("mqsi/EVAL_BSPLINE.f90")
spline = fmodpy.fimport("mqsi/SPLINE.f90")
mqsi = fmodpy.fimport("mqsi/MQSI.f90")
# Define an object that has all of the expected functions.
class splines:
    mqsi = mqsi.mqsi
    fit_spline = spline.fit_spline
    eval_spline = spline.eval_spline
    eval_bspline = bspline.eval_bspline


f = lambda x: x**2
interval = [0, 1]
size = 100

import numpy as np
x = np.linspace(*interval,size)
y = np.array([[f(v) for i in range(len(f.values[0]))]
              for v in x], dtype=float, order='F')
nb = len(y)
ncc = len(y[0])
print('x: ', x.shape)
print('y: ', y.shape)

x = np.array(x, dtype=float)
y = np.array(y, dtype=float, order='F')
sk = np.ones(nb*ncc + 2*ncc, dtype=float)
sc = np.ones(nb*ncc, dtype=float)

from util.system import Timer
t = Timer(); t.start()
sk, sc, info = splines.fit_spline(x, y, sk, sc)
t.stop()
print("t(): ",t())
print()
print("sk:   ", list(sk))
print("sc:   ", list(sc))
print("info: ", info)
assert( (info == 0) or (info > 10) )


