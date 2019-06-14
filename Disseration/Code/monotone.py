# Need to figure out an appropriate method for obtaianing a quintic
# spline. If an adjustment to one segment causes a neighbor to break,
# then we have to shrink the derivatives and repeat. There should be
# an exact value to which we must shrink the derivatives for
# monotonicity to be achievable. I'll have to algebraically solve for
# that.

# I need to identify the rescaling for DDX by algebraically solving
# given fixed DDX vector


from polynomial import fit

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

# Given x and y values, construct a monotonic quintic spline fit.
def monotone_cubic_spline(x, y):
    from polynomial import Spline
    f = fit(x, y, continuity=1)
    values = f.values
    # Make all pieces monotone.
    for i in range(len(values)-1):
        monotone_cubic_piece([x[i], x[i+1]], [values[i], values[i+1]])
    # Construct a new spline over the (updated) values for derivatives.
    return Spline(f.knots, f.values)

# Given a (x1, x2) and ([y1, d1y1], [y2, d1y2]), compute the rescaled
# y values for a monotone cubic piece.
def monotone_cubic_piece(x,y):
    # Compute the secant slope, the left slope ratio and the
    # right slope ratio for this interval of the function.
    secant_slope = (y[1][0] - y[0][0]) / (x[1] - x[0])
    left_ratio = y[0][1] / secant_slope
    right_ratio = y[1][1] / secant_slope
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
    changed = False
    # Extract local variable names from the provided points and
    # derivative information (to match the source paper).
    U0, U1 = x
    X0, DX0, DDX0 = y[0]
    X1, DX1, DDX1 = y[1]
    print()
    print("-"*70)
    print("        ","   DX0      DX1     DDX0     DDX1 ")
    print("Initial:",("%7.2f  "*4)%(DX0, DX1, DDX0, DDX1))
    v = X1 - X0
    h = (U1 - U0) / 2
    # Functions for evaluating constants in expressions (involving values)
    two_h_over_v = 2 * h / v
    h_sq_over_v = h * (h / v)
    A = lambda: DX0 * two_h_over_v
    B = lambda: DX1 * two_h_over_v
    C = lambda: DDX0 * h_sq_over_v 
    D = lambda: DDX1 * h_sq_over_v
    # Set DX0 and DX1 to be the median of these three choices.
    DX0 = sorted([0, DX0, 7*v / h])[1]
    DX1 = sorted([0, DX1, 7*v / h])[1]
    # Set A = max(0, A) and B = max(0, B).
    if A() < 0:
        DX0 = 0
        changed = True
    if B() < 0:
        DX1 = 0
        changed = True
    # Use a monotone cubic over this region if AB = 0.
    if (A()*B() == 0): return monotone_cubic_piece(x, y)
    # Scale the derivative vector to make tau1 positive.
    def tau_1(): return 24 + 2*(A()*B())**(1/2) - 3*(A()+B())
    if (tau_1() <= 0):
        a = A()
        b = B()
        # Compute the rescale factor necessary to make tau_1 0.
        rescale = 24 * (3*(a+b) + 2*(a*b)**(1/2)) / (9*(a**2+b**2) + 14*(a*b))
        # Shrink that rescale factor to ensure tau_1 becomes positive.
        sqrt_eps = 1.4901161193847656*10**(-8) 
        #          ^^ SQRT(EPSILON( 0.0_REAL64 ))
        rescale -= sqrt_eps**(1/2)
        # Rescale the derivatives
        DX0 *= rescale
        DX1 *= rescale
        changed = True
    assert(tau_1() > 0)
    # Compute DDX0 and DDX1 that satisfy monotonicity.


    # Scale (C,D) down until montonicity is achieved by using a simple verification.
    def is_monotone():
        a = B()
        b = 4*(B() - D())
        c = 30 - 6*(D() + C() + 2*B() + 2*A())
        d = 4*(A() + C())
        e = A()
        alpha = b * a**(-3/4) * e**(-1/4)
        beta  = c * a**(-1/2) * e**(-1/2)
        gamma = d * a**(-1/4) * e**(-3/4)
        if beta <= 6: small = - (beta + 2) / 2
        else:         small = -2 * (beta - 2)**(1/2)
        return (alpha > small) and (gamma > small)
    # Move the second derivative towards a working value until
    # monotonicity is achieved.
    target_DDX0 = (-(A()**(1/2) / 4) * (7*A()**(1/2) + 3*B()**(1/2)) * v) / h**2
    target_DDX1 = ( (B()**(1/2) / 4) * (3*A()**(1/2) + 7*B()**(1/2)) * v) / h**2
    steps = 0
    while (not is_monotone()):
        steps += 1
        DDX0 = .01*(target_DDX0) + .99*(DDX0)
        DDX1 = .01*(target_DDX1) + .99*(DDX1)
        if steps >= 1000: break
    print("Final:  ",("%7.2f  "*4)%(DX0, DX1, DDX0, DDX1))
    print("target_DDX0: ",target_DDX0)
    print("target_DDX1: ",target_DDX1)
    # Update "y" and return the updated version.
    y[0][:] = X0, DX0, DDX0
    y[1][:] = X1, DX1, DDX1
    changed = changed or (steps > 0)
    return changed

# --------------------------------------------------------------------

# 0, 4 -> Good
# 1, 4 -> Bad (now good)

SEED = 0
NODES = 13

# Generate random data to test the monotonic fit function.
import numpy as np
np.random.seed(SEED)
nodes = NODES + 2
x = np.linspace(0, 1, nodes)
y = sorted(np.random.normal(size=(nodes,)))
x -= min(x)
x /= max(x)
y -= min(y)
y /= max(y)
y -= max(y)
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

f = monotone_cubic_spline(x,y)
p.add_func("monotone c1", f, [min(x), max(x)], group='1m')
p.add_func("monotone c1 d1", f.derivative(), [min(x), max(x)], 
           dash="dash", color=p.color(p.color_num,alpha=.5), group='1m')

p.add_func("continuity 2", fit(x,y,continuity=2), [min(x), max(x)], group='2')
p.add_func("c2 d1", fit(x,y,continuity=2).derivative(), [min(x), max(x)], 
           dash="dash", color=p.color(p.color_num,alpha=.5), group='2')

f = monotone_quintic_spline(x,y)
p.add_func("monotone c2", f, [min(x), max(x)], group='2m', color=p.color(6))
p.add_func("monotone c2d1", f.derivative(), [min(x), max(x)], 
           dash="dash", color=p.color(6,alpha=.5), group='2m')

p.show(file_name="bad-hermite_interpolant.html")
