
TESTS = {}

# --------------------------------------------------------------------
# Watson test.
import numpy as np
from polynomial import Spline, fit
y = np.array([0,1,1,1,0,20,19,18,17,0,0,3,0,1,6,16,16.1,1], dtype=float)
x = np.linspace(0,1,len(y))
TESTS["watson_test"] = Spline(x, y[:,None])

# --------------------------------------------------------------------
# Sin function.
TESTS["sin_function"] = np.sin

# --------------------------------------------------------------------
# Large tangent test.
TESTS["large_tangent"] = lambda x: -(1.0 + ( 1.0 / (x-0.99) ))

# --------------------------------------------------------------------
# Random monotone data test.
TESTS["random_monotone"] = lambda x: np.random.random(size=len(x))


