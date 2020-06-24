

# Import the splines Fortran code.
import fmodpy

spline = fmodpy.fimport("SPLINE.f90")
# Define an object that has all of the expected functions.
class splines:
    fit_spline = spline.fit_spline
    eval_spline = spline.eval_spline

from polynomial import Spline
from fraction import Fraction
def exact(l):
    try: return [exact(v) for v in l]
    except TypeError: pass
    return Fraction(l)

# f = Polynomial([Fraction(1),Fraction(0),Fraction(0),Fraction(0)])
# f = Polynomial([Fraction(1),Fraction(0)])
interval = [0,1]
f = Spline(interval, exact([[1,-1],[1,1]]))
g1 = f.integral(c=0)
g2 = g1.integral(c=0)
g3 = g2.integral(c=0)
print("f: ",f)
print("g1: ",g1)
print("g2: ",g2)
print("g3: ",g3)
print("1/6: ",1/6)
print()

import numpy as np
x = np.linspace(*interval,2)
y = np.array([[f.derivative(i)(v) for i in range(len(f.values[0]))]
              for v in x], dtype=float, order='F')
nb = len(y)
ncc = len(y[0])
print('x: ', x.shape)
print('y: ', y.shape)

x = np.array(x, dtype=float)
y = np.array(y, dtype=float, order='F')
sk = np.ones(nb*ncc + 2*ncc, dtype=float)
sc = np.ones(nb*ncc+1, dtype=float)
sk, sc, info = splines.fit_spline(x, y, sk, sc)

print("sk:   ", list(sk))
print("sc:   ", list(sc))
print("info: ", info)
assert( (info == 0) or (info > 10) )



from util.plot import Plot

p = Plot()
p.add_func("f", f, interval)
p.add_func("g1", g1, interval)
p.add_func("g2", g2, interval)
p.add_func("g3", g3, interval)

x = np.linspace(*interval, 100)

y, info = splines.eval_spline(sk, sc, x.copy(), d=0) ; print(info)
p.add("sf", x, y, color=0, dash='dot')

y, info = splines.eval_spline(sk, sc, x.copy(), d=-1) ; print(info)
p.add("sg1", x, y, color=1, dash='dot')

y, info = splines.eval_spline(sk, sc, x.copy(), d=-2) ; print(info)
p.add("sg2", x, y, color=2, dash='dot')

y, info = splines.eval_spline(sk, sc, x.copy(), d=-3) ; print(info)
p.add("sg3", x, y, color=3, dash='dot')

p.show()
