
TESTS = {}

# --------------------------------------------------------------------
# Watson test.
from polynomial import Spline
y = [[0,1],[1,0],[1,0],[1,0],[0,0], # Plateu with curvature only on ends.
     [20,-1],[19,-1],[18,-1],[17,-1], # Sudden linear segment
     [0,0],[0,0],[3,0],[0,0], # Sudden flat with one peak
     [1,3],[6,9], # Exponential growth
     [16,.1],[16.1,.1], # small linear growth
     [1,-15]]
x = list(range(len(y)))
# Convert to cover the unit interval.
x = [v/(len(y)-1) for v in x]
for i in range(len(y)): y[i][1] *= (len(y)-1)
# Store the test.
TESTS["piecewise_polynomial"] = Spline(x, y)

# --------------------------------------------------------------------
# Signal function.
from numpy import sin, pi
TESTS["signal"] = lambda x: sin(4 * (2*pi) * x) / (x**2 + .1)

# --------------------------------------------------------------------
# Large tangent test.
TESTS["large_tangent"] = lambda x: -(1.0 + ( 1.0 / (x-0.99) ))

# --------------------------------------------------------------------
# Random data test.
from numpy import random
TESTS["random"] = lambda x: random.random(size=len(x)) if hasattr(x,"__len__") else random.random()

# --------------------------------------------------------------------
# Random monotone data test.
from numpy import random
TESTS["random_monotone"] = lambda x: sorted(random.random(size=len(x))) if hasattr(x,"__len__") else random.random()


