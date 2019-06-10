


# Given a (x1, x2) and ([y1, d1y1], [y2, d1y2]), compute the rescaled
# y values for a monotone cubic piece.
def monotone_cubic_piece(x,y):
    # Compute the secant slope, the left slope ratio and the
    # right slope ratio for this interval of the function.
    secant_slope = (y[1] - y[0]) / (x[1] - x[0])
    left_ratio = deriv[0] / secant_slope
    right_ratio = deriv[1] / secant_slope
    # ----------------------------------------------------------------
    #    USE PROJECTION ONTO CUBE FOR CUBIC SPLINE CONSTRUCTION
    # Determine which line segment it will project onto.
    mult = 3 / max(left_ratio, right_ratio)
    # Perform the projection onto the line segment by
    # shortening the vector to make the max coordinate 3.
    left_ratio = min(left_ratio, mult*left_ratio)
    right_ratio = min(right_ratio, mult*right_ratio)
    # ----------------------------------------------------------------
    # Set the derivative values based on (projected) monotone slope ratios.
    start_y = (tuple(y[0]), tuple(y[1]))
    y[0][1] = left_ratio * secant_slope
    y[1][1] = right_ratio * secant_slope
    changed = (start_y != (tuple(y[0]), tuple(y[1])))
    return changed

# Given a (x1, x2) and ([y1, d1y1, d2y1], [y2, d1y2, d2y2]), compute
# the rescaled derivative values for a monotone quintic piece.
def monotone_quintic_piece(x, y):
    U0, U1 = x
    X0, X1 = y[0][0], y[1][0]
    DX0, DX1 = y[0][1], y[1][1]
    DDX0, DDX1 = y[0][2], y[1][2]
    v = X1 - X0
    h = (U1 - U0) / 2
    A = lambda: (2 * h * DX0) / v
    B = lambda: (2 * h * DX1) / v
    C = lambda: (h*h * DDX0) / v
    D = lambda: (h*h * DDX1) / v
    tau1 = lambda: 24 + 2*(A()*B())**(1/2) - 3*(A() + B())
    alpha = lambda eta: (4 * (B() - eta[1])) / (A()**(1/4) * B()**(3/4))
    beta = lambda eta: ( 6*(eta[1] - eta[0] - 2*B() - 2*A() + 5) /
                         (A()**(1/2) + B()**(1/2)) )
    gamma = lambda eta: (4 * (A() + eta[0])) / (A()**(3/4) * B()**(1/4))
    eta_star = lambda eta, rho: (rho*eta[0] + (1 - rho)*eta[0],
                                 rho*eta[1] + (1 - rho)*eta[1])
    # Function for checking if a quintic is monotone based on the
    # coefficients of the (monomial) derivative.
    def is_monotone():
        a = B()
        b = 4*(B() - D())
        c = 6*(D() - C() - 2*B() - 2*A() + 5)
        d = 4*(A() + C())
        e = A()
        alpha = b * a**(-3/4) * e**(-1/4)
        beta = c * a**(-1/2) * e**(-1/2)
        gamma = d * a**(-1/4) * e**(-3/4)
        if beta <= 6: small = - (beta + 2) / 2
        else:         small = -2 * (beta - 2)**(1/2)
        return (alpha > small) and (gamma > small)

    # Set DX0 and DX1 to be the median of these three choices.
    DX0 = sorted([0, DX0, 7*v / h])[1]
    DX1 = sorted([0, DX1, 7*v / h])[1]
    # Set A = max(0, A) and B = max(0, B)
    if A() < 0: DX0 = 0
    if B() < 0: DX1 = 0
    # Use a monotone cubic over this region if AB = 0.
    if (A()*B() == 0): return monotone_cubic_piece(x, [[X0,X1], [DX0, DX1]])
    # Scale the derivative vector to be in the box (0,6], (0,6].
    multiplier = 6 / max(DX0, DX1)
    DX0 = min(DX0, multiplier*DX0)
    DX1 = min(DX1, multiplier*DX1)
    # Scale (C,D) down until montonicity is achieved.
    eta0 = (-(A()**(1/2) / 4) * (7*A()**(1/2) + 3*B()**(1/2)),
             (B()**(1/2) / 4) * (3*A()**(1/2) + 7*B()**(1/2)) )
    target_DDX0 = (eta0[0] * v) / h**2
    target_DDX1 = (eta0[1] * v) / h**2
    # Move the second derivative towards a working value until
    # monotonicity is achieved.
    steps = 0
    while (not is_monotone()):
        steps += 1
        DDX0 = .01*(target_DDX0) + .99*(DDX0)
        DDX1 = .01*(target_DDX1) + .99*(DDX1)
        if steps >= 1000: break
    # Update "y" and return the updated version.
    y[0][:] = X0, DX0, DDX0
    y[1][:] = X1, DX1, DDX1
    changed = (steps > 0)
    return changed

# Given data points "x" and data values "y", construct a 
# monotone interpolating spline over the given points with
# specified level of continuity.
def fit(x, y, continuity=0):
    from polynomial import fill_derivative, Spline
    knots = list(x)
    values = [[v] for v in y]
    # Construct further derivatives and refine the approximation
    # ensuring monotonicity in the process.
    for i in range(1,continuity+1):
        deriv = fill_derivative(knots, [v[-1] for v in values])
        for v,d in zip(values,deriv): v.append(d)
    # Return the interpolating spline.
    return Spline(knots, values)

# Given x and y values, construct a monotonic quintic spline fit.
def monotone_quintic_spline(x, y):
    from polynomial import Spline
    f = fit(x, y, continuity=2)
    values = f.values
    # Make all pieces monotone.
    for i in range(len(values)-1):
        monotone_quintic_piece([x[i], x[i+1]], [values[i], values[i+1]])

    # # Continue looping over the polynomial pieces until all pieces are
    # # monotone because fixing one may accidentally break the previous.
    # changed = True
    # while changed:
    #     changed = False
    #     for i in range(len(values)-1):
    #         changed = changed or monotone_quintic_piece(
    #             [x[i], x[i+1]], [values[i], values[i+1]])

    # Construct a new spline over the (updated) values for derivatives.
    return Spline(f.knots, f.values)

SEED = 3
NODES = 4

# Generate random data to test the monotonic fit function.
import numpy as np
np.random.seed(SEED)
nodes = NODES + 2
x = np.linspace(0, 1, nodes)
y = sorted(np.random.normal(size=(nodes,)))
# Convert these arrays to lists.
x, y = list(x), list(y)


# Generate a plot to see what it all looks like.
from util.plot import Plot
p = Plot()
p.add("Points", x, y)
p.add_func("continuity 0", fit(x,y,continuity=0), [min(x), max(x)], group='0')
p.add_func("c0 d1", fit(x,y,continuity=0).derivative(), [min(x), max(x)], 
           dash="dash", color=p.color(p.color_num,alpha=.5), group='0')
p.add_func("continuity 1", fit(x,y,continuity=1), [min(x), max(x)], group='1')
p.add_func("c1 d1", fit(x,y,continuity=1).derivative(), [min(x), max(x)], 
           dash="dash", color=p.color(p.color_num,alpha=.5), group='1')
p.add_func("continuity 2", fit(x,y,continuity=2), [min(x), max(x)], group='2')
p.add_func("c2 d1", fit(x,y,continuity=2).derivative(), [min(x), max(x)], 
           dash="dash", color=p.color(p.color_num,alpha=.5), group='2')
f = monotone_quintic_spline(x,y)
p.add_func("monotone c2", f, [min(x), max(x)], group='3')
p.add_func("monotone c2d1", f.derivative(), [min(x), max(x)], 
           dash="dash", color=p.color(p.color_num,alpha=.5), group='3')

p.show(file_name="hermite_interpolant.html")
