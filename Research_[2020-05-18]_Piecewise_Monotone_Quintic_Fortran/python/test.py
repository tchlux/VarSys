import numpy as np
# x =  [1., 2., 2.5, 3., 4.]
# y = [1., 4., 6.25, 7., 10.]
# y = [1., 4., 6.25, 6.999999999999999, 9.999999999999993]

# x = np.linspace(0,1,20)
# y = np.sin(8*np.pi*x) / (x**2 + 0.1)

np.random.seed(2)
N = 16
x = np.random.random(size=(N,))
y = np.random.random(size=(N,))
x.sort(); y.sort()

from mqsi import monotone_quintic_spline, weighted_harmonic, min_curve

q1 = monotone_quintic_spline(x,y)
q2 = monotone_quintic_spline(x,y, estimator=weighted_harmonic)
f = min_curve(x, y, t=10)

from util.plot import Plot
p = Plot()
p.add("Data", x, y)
p.add_func("Q", f, [min(x), max(x)])
p.add_func("QF", q1, [min(x), max(x)])
p.add_func("WH", q2, [min(x), max(x)])
p.show()



