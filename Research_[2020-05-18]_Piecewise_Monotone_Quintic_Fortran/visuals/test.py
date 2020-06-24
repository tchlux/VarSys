from util.plot import Plot
from mqsi import Polynomial, monotone_quintic_spline

p = Plot()

f1 = lambda x: x**2
f2 = lambda x, m=1: m*(x-2)**2 + 6

def linspace(l,r,n): return [l+(i/(n-1))*(r-l) for i in range(n)]

eps = 1 - 2**(-51)
x = linspace(1, 4, 5)
x[1] = 2*x[1]/3 + x[2]/3
x[3] = 2*x[3]/3 + x[2]/3
y1 = [f1(x[0]), f1(x[1]), f1(x[2]), f2(x[3]), f2(x[4])]
gap = f2(x[2]) - f2(x[2],eps)
y2 = [f1(x[0]), f1(x[1]), f1(x[2]), f2(x[3],eps)+gap, f2(x[4],eps)+gap]

q1 = monotone_quintic_spline(x, y1, monotone=True)
q2 = monotone_quintic_spline(x, y2, monotone=True)
# q1o = monotone_quintic_spline(x, y1, monotone=False)
# q2o = monotone_quintic_spline(x, y2, monotone=False)

p.add("data 1", x, y1, color=0)
p.add("data 2", x, y2, color=0)
bounds = [min(x), max(x)]
p.add_function("f1", f1, bounds, color=3, opacity=.2, dash="dash")
p.add_function("f2", f2, bounds, color=3, opacity=.2, dash="dash")
p.add_function("Q1", q1, bounds, color=1)
p.add_function("Q2", q2, bounds, color=2)
# p.add_function("Q1 qfm", q1o, bounds, color=1, opacity=.5)
# p.add_function("Q2 qfm", q2o, bounds, color=2, opacity=.5)
p.show()
