# This file provides the function `monotone_quintic_spline` which
# constructs spline interpolant over provided data. It starts by using
# the Newton method to construct approximations for the first and
# second derivative at each of the data points. Then it iteratively
# adjusts the derivatives and second derivatives until all quintic
# polynomial pieces are monotone. A `Spline` object is returned.

from polynomial import fit

# Given x and y values, construct a monotonic quintic spline fit.
def monotone_quintic_spline(x, y, ends=2, mids=2, fix_previous=True):
    from polynomial import Spline
    f = fit(x, y, continuity=2, non_decreasing=True, ends=ends, mids=mids)
    values = f.values
    # Make all pieces monotone.
    i = 0
    print()
    while i < len(values)-1:
        changed = nondecreasing_quintic(
            [x[i], x[i+1]], [values[i], values[i+1]])
        if fix_previous and changed and (i > 0):
            step = 0
            # While the previous interval is also broken,
            while not is_quintic_nondecreasing(
                    [x[i-1], x[i]], [values[i-1], values[i]]):
                step += 1
                # Shrink the derivative value at this point.
                values[i][1] *= .99

                # Fix the previous interval.
                nondecreasing_quintic(
                    [x[i-1], x[i]], [values[i-1], values[i]])
                # Fix this interval again.
                nondecreasing_quintic(
                    [x[i], x[i+1]], [values[i], values[i+1]])
                if step >= 100:
                    print()
                    print('-'*70)
                    print()
                    print(f"Interval {i}: {x[i]:.3f} {x[i+1]:.3f}")
                    print()
                    print('-'*70)
                    print()
                    raise(Exception("Error for this interval.."))
        # Increment to work on the next interval.
        i += 1
    # Construct a new spline over the (updated) values for derivatives.
    return Spline(f.knots, f.values)

# Given x and y values, construct a monotonic quintic spline fit.
def monotone_cubic_spline(x, y):
    from polynomial import Spline
    f = fit(x, y, continuity=1, non_decreasing=True)
    values = f.values
    # Make all pieces monotone.
    for i in range(len(values)-1):
        nondecreasing_cubic([x[i], x[i+1]], [values[i], values[i+1]])
    # Construct a new spline over the (updated) values for derivatives.
    return Spline(f.knots, f.values)

# Given a (x1, x2) and ([y1, d1y1], [y2, d1y2]), compute the rescaled
# y values for a monotone cubic piece.
def nondecreasing_cubic(x,y):
    # Compute the secant slope, the left slope ratio and the
    # right slope ratio for this interval of the function.
    secant_slope = (y[1][0] - y[0][0]) / (x[1] - x[0])
    A = y[0][1] / secant_slope # (left slope ratio)
    B = y[1][1] / secant_slope # (right slope ratio)
    # ----------------------------------------------------------------
    #    USE PROJECTION ONTO CUBE FOR CUBIC SPLINE CONSTRUCTION
    # Determine which line segment it will project onto.
    mult = 3 / max(A, B)
    if (mult < 1):
        # Perform the projection onto the line segment by
        # shortening the vector to make the max coordinate 3.
        y[0][1] *= mult
        y[1][1] *= mult
        return True
    # ----------------------------------------------------------------
    return False

# Given a (x1, x2) and ([y1, d1y1, d2y1], [y2, d1y2, d2y2]), compute
# the rescaled derivative values for a monotone quintic piece.
def nondecreasing_quintic(x, y):
    changed = False
    # Extract local variable names from the provided points and
    # derivative information (to match the source paper).
    U0, U1 = x
    X0, DX0, DDX0 = y[0]
    X1, DX1, DDX1 = y[1]
    function_change = X1 - X0
    interval_width = U1 - U0
    interval_slope = function_change / interval_width
    # Set DX0 and DX1 to be the median of these three choices.
    DX0 = sorted([0, DX0, 14*interval_slope])[1]
    DX1 = sorted([0, DX1, 14*interval_slope])[1]
    # Compute repeatedly used values "A" (left ratio) and "B" (right ratio).
    A = DX0 / interval_slope
    B = DX1 / interval_slope
    # Set A = max(0, A) and B = max(0, B).
    if (A < 2**(-26)): DX0 = A = 0
    if (B < 2**(-26)): DX1 = B = 0
    # Use a monotone cubic over this region if AB = 0.
    if (A*B < 2**(-26)):
        y[0][2] = y[1][2] = 0
        return True or nondecreasing_cubic(x, y)
    # Scale the derivative vector to make tau_1 positive.
    tau_1 = 24 + 2*(A*B)**(1/2) - 3*(A+B)
    if (tau_1 <= 0):
        # Compute the rescale factor necessary to make tau_1 0.
        rescale = 24 * (3*(A+B) + 2*(A*B)**(1/2)) / (9*(A**2+B**2) + 14*(A*B))
        # Shrink that rescale factor to ensure tau_1 becomes positive.
        rescale -= rescale * 2**(-26) # SQRT(EPSILON( 0.0_REAL64 ))
        # Rescale the derivatives
        DX0 *= rescale
        DX1 *= rescale
        # Recompute A and B.
        A = DX0 / interval_slope
        B = DX1 / interval_slope
        # Record the change.
        changed = True
    # Make sure that the condition is met.
    tau_1 = 24 + 2*(A*B)**(1/2) - 3*(A+B)
    assert(tau_1 > 0)

    # Clip values that are too large (to ensure that the neighboring
    # regions can be adjusted without worry).
    mult = 8 / max(A, B)
    if (mult < 1):
        DX0 *= mult
        DX1 *= mult
        A = DX0 / interval_slope
        B = DX1 / interval_slope
        changed = True

    # Compute DDX0 and DDX1 that satisfy monotonicity by scaling (C,D)
    # down until montonicity is achieved (using binary search).
    alpha_constant = 4 * (B**(1/4) / A**(1/4))
    beta_constant  = (12 * (DX0+DX1) * (U1-U0) + 30 * (X0-X1)) / ((X0-X1) * A**(1/2) * B**(1/2))
    gamma_constant = 4 * ((DX0 * B**(3/4)) / (DX1 * A**(3/4)))
    alpha_multiplier = ((U0-U1) / DX1) * B**(1/4) / A**(1/4)
    beta_multiplier  = (3 * (U0-U1)**2) / (2 * (X0-X1) * A**(1/2) * B**(1/2))
    gamma_multiplier = ((U1-U0) / DX1) * B**(3/4) / A**(3/4)
    def is_monotone():
        a = alpha_constant + alpha_multiplier * DDX1
        g = gamma_constant + gamma_multiplier * DDX0
        b = beta_constant  + beta_multiplier  * (DDX0 - DDX1)
        if b <= 6: bound = - (b + 2) / 2
        else:      bound = -2 * (b - 2)**(1/2)
        return (a > bound) and (g > bound)

    # Move the second derivative towards a working value until
    # monotonicity is achieved.
    target_DDX0 = - A**(1/2) * (7*A**(1/2) + 3*B**(1/2)) * interval_slope / interval_width
    target_DDX1 =   B**(1/2) * (3*A**(1/2) + 7*B**(1/2)) * interval_slope / interval_width

    # If this function is not monotone, perform a binary
    # search for the nearest-to-original monotone DDX values.
    if not is_monotone():
        accuracy = 2**(-26) # = SQRT(EPSILON( 0.0_REAL64 ))
        original_DDX0, original_DDX1 = DDX0, DDX1
        low_bound,     upp_bound     = 0,    1
        # Continue dividing the interval in 2 to find the smallest
        # "upper bound" (amount target) that satisfies monotonicity.
        while ((upp_bound - low_bound) > accuracy):
            to_target = (upp_bound + low_bound) / 2
            # If we found the limit of floating point numbers, break.
            if ((to_target == upp_bound) or (to_target == low_bound)): break
            # Recompute DDX0 and DDX1 based on the to_target.
            DDX0 = (1-to_target) * original_DDX0 + to_target * target_DDX0
            DDX1 = (1-to_target) * original_DDX1 + to_target * target_DDX1
            # Otherwise, proceed with a binary seaarch.
            if is_monotone(): upp_bound = to_target
            else:             low_bound = to_target
        # Store the smallest amount "target" for DDX0 and DDX1 that is nondecreasing.
        DDX0 = (1-upp_bound) * original_DDX0 + upp_bound * target_DDX0
        DDX1 = (1-upp_bound) * original_DDX1 + upp_bound * target_DDX1
        changed = True
        # Verify that the function is monotone (according to check function).
        try: assert(is_monotone())
        except:
            print("function_change: ",function_change)
            print("interval_width:  ",interval_width)
            print("interval_slope:  ",interval_slope)
            print("x: ",x)
            print("y: ",y)
            print("DX0: ",DX0)
            print("DX1: ",DX1)
            print("DDX0: ",DDX0)
            print("DDX1: ",DDX1)
            print("target_DDX0: ",target_DDX0)
            print("target_DDX1: ",target_DDX1)
            print("low_bound: ",low_bound)
            print("upp_bound: ",upp_bound)
            print()
            a = alpha_constant + alpha_multiplier * DDX1
            g = gamma_constant + gamma_multiplier * DDX0
            b = beta_constant  + beta_multiplier  * (DDX0 - DDX1)
            print("a: ",a)
            print("g: ",g)
            print("b: ",b)
            print()
            raise(Exception("Monotonicity was violated for unknown reasons."))

    # Update "y" and return the updated version.
    y[0][:] = X0, DX0, DDX0
    y[1][:] = X1, DX1, DDX1
    return changed

# Function for testing if a quintic piece is monotone.
def is_quintic_nondecreasing(x, y):
    # Extract meaningful components of provided variables.
    U0, U1 = map(float,x)
    X0, DX0, DDX0 = map(float,y[0])
    X1, DX1, DDX1 = map(float,y[1])
    # Compute useful intermediate values.
    function_change = X1 - X0
    interval_width = U1 - U0
    interval_slope = function_change / interval_width
    # Compute left and right slope ratios.
    A = DX0 / interval_slope
    B = DX1 / interval_slope
    if (A*B == 0): return True
    # Compute constants and multipliers for determining monotonicity.
    alpha_constant = 4 * (B**(1/4) / A**(1/4))
    beta_constant  = (12 * (DX0+DX1) * (U1-U0) + 30 * (X0-X1)) / ((X0-X1) * A**(1/2) * B**(1/2))
    gamma_constant = 4 * ((DX0 * B**(3/4)) / (DX1 * A**(3/4)))
    alpha_multiplier = ((U0-U1) / DX1) * B**(1/4) / A**(1/4)
    beta_multiplier  = (3 * (U0-U1)**2) / (2 * (X0-X1) * A**(1/2) * B**(1/2))
    gamma_multiplier = ((U1-U0) / DX1) * B**(3/4) / A**(3/4)
    # Compute monotonicity.
    a = alpha_constant + alpha_multiplier * DDX1
    g = gamma_constant + gamma_multiplier * DDX0
    b = beta_constant  + beta_multiplier  * (DDX0 - DDX1)
    if b <= 6: bound = - (b + 2) / 2
    else:      bound = -2 * (b - 2)**(1/2)
    return (a > bound) and (g > bound)


if __name__ == "__main__":
    # --------------------------------------------------------------------

    # 0, 4 -> Good
    # 1, 4 -> Bad (now good)
    # 0, 13 -> Bad (now good)
    # 0, 100 -> Bad

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

    p.add_func("continuity 0", fit(x,y,continuity=0, non_decreasing=True),
               [min(x), max(x)], group='0')
    p.add_func("c0 d1", fit(x,y,continuity=0, non_decreasing=True).derivative(),
               [min(x), max(x)], dash="dash", color=p.color(p.color_num,alpha=.5), group='0')

    p.add_func("continuity 1", fit(x,y,continuity=1,
                                   non_decreasing=False, ends=1, mids=0), 
               [min(x), max(x)], group='1')
    p.add_func("c1 d1", fit(x,y,continuity=1, non_decreasing=True).derivative(), 
               [min(x), max(x)], dash="dash", color=p.color(p.color_num,alpha=.5), group='1')

    f = monotone_cubic_spline(x,y)
    p.add_func("monotone c1", f, [min(x), max(x)], group='1m')
    p.add_func("monotone c1 d1", f.derivative(), [min(x), max(x)], 
               dash="dash", color=p.color(p.color_num,alpha=.5), group='1m')

    p.add_func("continuity 2", fit(x,y,continuity=2, non_decreasing=True), 
               [min(x), max(x)], group='2')
    p.add_func("c2 d1", fit(x,y,continuity=2, non_decreasing=True).derivative(), 
               [min(x), max(x)], dash="dash", color=p.color(p.color_num,alpha=.5), group='2')

    f = monotone_quintic_spline(x,y, fix_previous=False)
    p.add_func("monotone c2 (no fix)", f, [min(x), max(x)], group='2m', color=p.color(6))
    p.add_func("monotone c2d1 (no fix)", f.derivative(), [min(x), max(x)], 
               dash="dash", color=p.color(6,alpha=.5), group='2m')

    print("Making full working spline..")
    f = monotone_quintic_spline(x,y)
    print("done.")
    p.add_func("monotone c2", f, [min(x), max(x)], group='2mf', color=p.color(7))
    p.add_func("monotone c2d1", f.derivative(), [min(x), max(x)], 
               dash="dash", color=p.color(7,alpha=.5), group='2mf')

    p.show(file_name="monotone_quintic_interpolating_spline.html")
