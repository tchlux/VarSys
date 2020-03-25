

# Configuration.
n = 4
lower = 0
upper = n


# Construct a point set.
from numpy import linspace
x = linspace(lower, upper, n+1)

# Construct a plot.
from util.plot import Plot
p = Plot()
# Construct the lagrange polynomials.
for k in range(n+1):
    from util.math import product
    def lagrange_polynomial(z):
        return product( (z - x[i]) / (x[k] - x[i]) for i in range(n+1) if i != k)
    # add to the plot.
    p.add_func(f"Lagrange polynomial {k}", lagrange_polynomial, [lower, upper])
p.show()

