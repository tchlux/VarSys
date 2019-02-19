# Differentiable vector operators.
from autograd import numpy as np
from autograd import elementwise_grad as deriv
from autograd.numpy.linalg import norm
# My utilities.
from util.algorithms import Delaunay
from util.system import save, load
from util.plot import Plot
from util.data import Data

# Seed the random number generator.
seed = 1
np.random.seed(seed)
tests = 6000
d = 1
n = 50
low = 0
upp = 1

# low = 0.0000000000001
# upp = 1 / norm(np.ones(d))

# Define a tricky sin function.
# divisor = d * norm(np.ones(d))
# f = lambda x: np.sum(np.sin(x)) / divisor

# mult =  (1*np.pi)
# f = lambda x: np.sin(norm(x) * mult) / mult**2
# g = deriv(f)
# h = lambda x: norm(g(x))

def f(x):
    x = x * (4 / upp)
    out = np.where( x <= 4, (x**2)/2 - 4*x + 8, 1 )
    out = np.where( x <= 3, -(x**2)/2 + 2*x - 1, out )
    out = np.where( x <= 1, (x**2)/2, out )
    out = np.where( x <= 0, 0, out )
    return np.sum(out) / (4/upp)**2

g = deriv(f)

# v = 0
# print(f(np.ones(d)*v))
# exit()

# print()
# print(f"{f(np.ones(d)*.5):.2f} {f(np.ones(d)*low):.2f} {f(np.ones(d)*upp):.2f}")
# if (d == 2):
#     print(f"[%.2f %.2f] [%.2f %.2f] [%.2f %.2f]"%(
#         tuple(g(np.ones(d)*.5)) + tuple(g(np.ones(d)*low)) + tuple(g(np.ones(d)*upp))))
# print(f"{h(np.ones(d)*.5):.2f} {h(np.ones(d)*low):.2f} {h(np.ones(d)*upp):.2f}")
# print()

# Define the error bound functions.
lipschitz_bound = lambda x: x
diagonal = np.ones(d)*(1/(n-1))
theorem_bound = lambda x: (x**2)/2 + d**(1/2) * ((upp - low) / (n-1)) * x

# Generate the grid of evaluation points.
train_x = np.vstack(x.flatten() for x in np.meshgrid(*(np.linspace(low,upp,n),)*d)).T
# Evaluate at all training pointss.
train_y = np.array([f(x) for x in train_x])
# Construct and fit a Delaunay model.
model = Delaunay(parallel=True)
model.fit(train_x, train_y)

print("Evaluating model at test points..")
try: test_x, test_y, errors, dists = load(None)
except:
    # Construct test points and measure error.
    test_x = (np.random.random(size=(tests,d)) * (upp - low)) - low
    test_y = np.array([f(x) for x in test_x])
    errors = np.array([abs(float(v)) for v in (model(test_x) - test_y)])
    dists =  np.array([np.min(norm(x-train_x,axis=1)) for x in test_x])
    save((test_x, test_y, errors, dists))


print("Making plot of results..")
p = Plot("", "Distance to Nearest Vertex", "Absolute Error")
p.add("Errors", dists, errors)
p.add_func("Lipschitz Error Bound", lipschitz_bound, [0,max(dists)])
x = .035 * (max(dists) - min(dists)) + min(dists)
p.add_annotation("Lipschitz Error Bound", x, lipschitz_bound(x), ax=(x/2))
p.add_func("Theorem Error Bound", theorem_bound, [0,max(dists)])
x = .07 * (max(dists) - min(dists)) + min(dists)
p.add_annotation("Theorem Error Bound", x, theorem_bound(x), ax=(2/3)*x)
p.show(y_range=[0,max(errors)], show_legend=False, show=(d>2))

if (d <= 2):
    p = Plot()
    p.add("Interpolation Points", *(train_x.T), train_y)
    p.add_func("Function", f, *([low,upp],)*d)
    p.add_func("Derivative", g, *([low,upp],)*d)
    p.show(append=True)
