# - Setting first and second derivatives to zero causes all subsequent
#   intervals to be difficult to fix. Look into why this is happening,
#   involves the simplified conditions.


# - Identify situation in which fixing one interval breaks previous.
# - Identify if, post-fix, an interval can be the left side of ^^
# - Guarantee that fixing an interval cannot break interval to right.

# - identify a scaling vector that can be applied to the derivative
#   and second derivative at one end of an interval that is guaranteed
#   not to break monotonicity.
# - prove that conditions on two neighboring intervals can arrive at
#   contradiction changes in the value of the first derivative and /
#   or second derivative.


# This file provides the function `monotone_quintic_spline` which
# constructs spline interpolant over provided data. It starts by using
# the Newton method to construct approximations for the first and
# second derivative at each of the data points. Then it iteratively
# adjusts the derivatives and second derivatives until all quintic
# polynomial pieces are monotone. A `Spline` object is returned.

from polynomial import fit

# Given x and y values, construct a monotonic quintic spline fit.
def monotone_quintic_spline(x, y, ends=2, mids=2, fix_previous=True, exact=True):
    print()
    from polynomial import Spline
    if exact:
        from fraction import Fraction
        x = list(map(Fraction, x))
        y = list(map(Fraction, y))
    # Construct an initial fit that is twice continuous.
    f = fit(x, y, continuity=2, non_decreasing=True, ends=ends, mids=mids)
    values = f.values
    # Make all pieces monotone.
    i = 0
    while i < len(values)-1:
        changed = monotone_quintic(
            [x[i], x[i+1]], [values[i], values[i+1]])
        # Maintain exactness if necessary!
        if exact:
            values[i][:] = map(Fraction, values[i])
            values[i+1][:] = map(Fraction, values[i+1])

        if fix_previous and changed and (i > 0):
            step = 0
            changed = True

            # While the previous interval is also broken,
            while changed:
                step += 1
                # Shrink the derivative value at the point between intervals.
                values[i][1] *= 9
                values[i][1] /= 10
                # Shrink the derivative value on the right side of this interval.
                values[i+1][1] *= 9
                values[i+1][1] /= 10
                # Fix the previous interval.
                changed = monotone_quintic(
                    [x[i-1], x[i]], [values[i-1], values[i]])
                if exact:
                    values[i-1][:] = map(Fraction, values[i-1])
                    values[i][:] = map(Fraction, values[i])
                # Fix this interval again.
                changed = monotone_quintic(
                    [x[i], x[i+1]], [values[i], values[i+1]]) or changed
                if exact:
                    values[i][:] = map(Fraction, values[i])
                    values[i+1][:] = map(Fraction, values[i+1])
                if step >= 500:
                    # Convert fractions into floats for formatted printing.
                    xl = float(x[i-1])
                    xm = float(x[i])
                    xr = float(x[i+1])
                    vl = list(map(float,values[i-1]))
                    vm = list(map(float,values[i]))
                    vr = list(map(float,values[i+1]))
                    print()
                    print('-'*70)
                    print()
                    print(f"Bad interval {i+1}, points {i+1} to {i+2}: ({xm:.3f}, {xr:.3f})")
                    print(f"  {xl}: {vl}")
                    print(f"  {xm}: {vm}")
                    print(f"  {xr}: {vr}")
                    print()
                    print('-'*70)
                    print()
                    values[i][1] *= 0
                    values[i][2] *= 0
                    values[i+1][1] *= 0
                    values[i+1][2] *= 0
                    break
                    # raise(Exception("Error for this interval.."))
            print(f"Interval {i+1} changed ({step} corrections)")
        elif (not changed):
            print(f"Interval {i+1} ...")
        else:
            print(f"Interval {i+1} changed")
        # Increment to work on the next interval.
        i += 1
    print()
    # Construct a new spline over the (updated) values for derivatives.
    return Spline(f.knots, f.values)

# Given x and y values, construct a monotonic quintic spline fit.
def monotone_cubic_spline(x, y):
    from polynomial import Spline
    f = fit(x, y, continuity=1, non_decreasing=True)
    values = f.values
    # Make all pieces monotone.
    for i in range(len(values)-1):
        monotone_cubic([x[i], x[i+1]], [values[i], values[i+1]])
    # Construct a new spline over the (updated) values for derivatives.
    return Spline(f.knots, f.values)

# Given a (x1, x2) and ([y1, d1y1], [y2, d1y2]), compute
# the rescaled y values to make this piece monotone. 
def monotone_cubic(x, y):
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
# the rescaled y values to make this piece monotone. 
def monotone_quintic_simplified(x, y):
    # Use simplified conditions from 1988 positivity conditions for
    # quartic polynomials to enforce monotonicity on this interval.
    U0, U1 = x
    X0, DX0, DDX0 = y[0]
    X1, DX1, DDX1 = y[1]
    w = U1 - U0
    changed = False
    # Enforce \gamma >= \delta
    min_DDX0 = -6 * DX0 / w
    if (min_DDX0 > DDX0):
        DDX0 = min_DDX0
        y[0][2] = DDX0
        changed = True
    # Enforce \alpha >= 0
    min_DDX1 = (14*DX0 + 16*DX1) / w  +  DDX0  +  60 * (X1 - X0) / (w*w)
    if (min_DDX1 > DDX1):
        DDX1 = min_DDX1
        y[1][2] = DDX1
        changed = True
    # Enforce \beta >= \alpha
    min_DDX1 = 6 * DX0 - 4 * DX1 - DDX0
    if (min_DDX1 > DDX1):
        y[1][2] = min_DDX1
        changed = True
    # Return whether or not the values were changed.
    return changed


# Given a (x1, x2) and ([y1, d1y1, d2y1], [y2, d1y2, d2y2]), compute
# the rescaled derivative values for a monotone quintic piece.
def monotone_quintic(x, y):
    changed = False
    # Extract local variable names from the provided points and
    # derivative information (to match the source paper).
    U0, U1 = x
    X0, DX0, DDX0 = y[0]
    X1, DX1, DDX1 = y[1]
    function_change = X1 - X0
    interval_width = U1 - U0
    interval_slope = function_change / interval_width
    # Handle an unchanging interval.
    if (interval_slope == 0):
        changed = any(v != 0 for v in (DX0, DX1, DDX0, DDX1))
        y[0][:] = X0, 0, 0
        y[1][:] = X1, 0, 0
        return changed
    sign = (-1) ** int(function_change < 0)
    # Set DX0 and DX1 to be the median of these three choices.
    if not (0 <= sign*DX0 <= sign*14*interval_slope):
        new_val = sign * (sorted([0, sign*DX0, sign*14*interval_slope])[1])
        changed = (new_val != DX0)
        if changed: DX0 = new_val
    if not (0 <= sign*DX1 <= sign*14*interval_slope):
        new_val = sign * (sorted([0, sign*DX1, sign*14*interval_slope])[1])
        changed = (new_val != DX1)
        if changed: DX1 = new_val
    # Compute repeatedly used values "A" (left ratio) and "B" (right ratio).
    A = DX0 / interval_slope
    B = DX1 / interval_slope
    assert(A >= 0)
    assert(B >= 0)
    # Use a monotone cubic over this region if AB = 0.
    if (A*B <= 0): return monotone_quintic_simplified(x, y) or True

    # Clip derivative values that are too large (to ensure that
    # shrinking the derivative vectors on either end will not break
    # monotonicity). (clipping at 6 alone is enough, with 8 needs more)
    mult = (6 / max(A, B))
    if (mult < 1):
        DX0 *= mult
        DX1 *= mult
        A = DX0 / interval_slope
        B = DX1 / interval_slope
        changed = True
    # Make sure that the first monotonicity condition is met.
    tau_1 = 24 + 2*(A*B)**(1/2) - 3*(A+B)
    try: assert(tau_1 >= 0)
    except:
        class NonMonotoneTau(Exception): pass
        raise(NonMonotoneTau(f"Bad Tau 1 value: {tau_1}\n A: {A}\n B: {B}"))

    # Compute DDX0 and DDX1 that satisfy monotonicity by scaling (C,D)
    # down until montonicity is achieved (using binary search).
    alpha_constant   = 4 * (B**(1/4) / A**(1/4))
    alpha_multiplier = ((U0-U1) / DX1) * B**(1/4) / A**(1/4)
    gamma_constant   = 4 * (DX0 / DX1) * (B**(3/4) / A**(3/4))
    gamma_multiplier = ((U1-U0) / DX1) * (B**(3/4) / A**(3/4))
    beta_constant    = (12 * (DX0+DX1) * (U1-U0) + 30 * (X0-X1)) / ((X0-X1) * A**(1/2) * B**(1/2))
    beta_multiplier  = (3 * (U0-U1)**2) / (2 * (X0-X1) * A**(1/2) * B**(1/2))
    def is_monotone():
        a = alpha_constant + alpha_multiplier * DDX1
        g = gamma_constant + gamma_multiplier * DDX0
        b = beta_constant  + beta_multiplier  * (DDX0 - DDX1)
        if b <= 6: bound = - (b + 2) / 2
        else:      bound = -2 * (b - 2)**(1/2)
        return (a+2**(-26) > bound) and (g+2**(-26) > bound)

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
            to_target = upp_bound / 2  +  low_bound / 2
            # If we found the limit of floating point numbers, break.
            if ((to_target == upp_bound) or (to_target == low_bound)): break
            # Recompute DDX0 and DDX1 based on the to_target.
            DDX0 = (1-to_target) * original_DDX0 + to_target * target_DDX0
            DDX1 = (1-to_target) * original_DDX1 + to_target * target_DDX1
            # Otherwise, proceed with a binary seaarch.
            if is_monotone(): upp_bound = to_target
            else:             low_bound = to_target
        # Store the smallest amount "target" for DDX0 and DDX1 that is monotone.
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
            if b <= 6: bound = - (b + 2) / 2
            else:      bound = -2 * (b - 2)**(1/2)
            print("a:     ",a)
            print("g:     ",g)
            print("b:     ",b)
            print("bound: ",bound)
            print()
            raise(Exception("Monotonicity was violated for unknown reasons."))

    # Update "y" and return the updated version.
    y[0][:] = X0, DX0, DDX0
    y[1][:] = X1, DX1, DDX1
    return changed


if __name__ == "__main__":
    # --------------------------------------------------------------------
    #               TEST CASES
    # 
    # 0, 4, (None)      -> Good
    # 1, 4, (None)      -> Bad (not monotone) (now good) [bad first derivative conditions]
    # 0, 13, (None)     -> Bad (not monotone after first pass) (now good) [no previous-interval fix]
    # 0, 30, (-3,None)  -> Bad (far right still not monotone after passes) (now good) [bad cubic usage]
    # 0, 100, (-12,-7)  -> Bad (fixing previous breaks second previous) (now good) [bad while loop condition]
    # 0, 100, (-4,None) -> Bad (last interval is not monotone) (now good) [forced first and second derivative to 0 when failing]
    # 0, 1000, (None)   -> Ehh (once one interval zeros, rest zero)
    # 
    # --------------------------------------------------------------------

    # SEED = 0
    # NODES = 100
    # SUBSET = slice(None) 

    SEED = 0
    NODES = 4
    SUBSET = slice(None)

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
    x, y = list(x)[SUBSET], list(y)[SUBSET]

    # Convert these arrays to exact arithmetic.
    from util.math import Fraction
    x = list(map(Fraction, x))
    y = list(map(Fraction, y))
    interval = [float(min(x)), float(max(x))]

    # Generate a plot to see what it all looks like.
    from util.plot import Plot
    p = Plot()
    p.add("Points", list(map(float, x)), list(map(float, y)))

    kwargs = dict(plot_points=1000)
    # Continuity 0
    p.add_func("continuity 0", fit(x,y,continuity=0, non_decreasing=True),
               interval, group='0', **kwargs)
    p.add_func("c0 d1", fit(x,y,continuity=0, non_decreasing=True).derivative(),
               interval, dash="dash", color=p.color(p.color_num,alpha=.5), 
               group='0', **kwargs)
    # Continuity 1
    p.add_func("continuity 1", fit(x,y,continuity=1,
                                   non_decreasing=False, ends=1, mids=0), 
               interval, group='1', **kwargs)
    p.add_func("c1 d1", fit(x,y,continuity=1, non_decreasing=True).derivative(), 
               interval, dash="dash", color=p.color(p.color_num,alpha=.5), 
               group='1', **kwargs)
    # Continuity 2
    f = monotone_cubic_spline(x,y)
    p.add_func("monotone c1", f, interval, group='1m', **kwargs)
    p.add_func("monotone c1 d1", f.derivative(), interval, 
               dash="dash", color=p.color(p.color_num,alpha=.5), group='1m', **kwargs)

    p.add_func("continuity 2", fit(x,y,continuity=2, non_decreasing=True), 
               interval, group='2', **kwargs)
    p.add_func("c2 d1", fit(x,y,continuity=2, non_decreasing=True).derivative(), 
               interval, dash="dash", color=p.color(p.color_num,alpha=.5), group='2', **kwargs)
    # Continuity with one-pass monotone fix.
    f = monotone_quintic_spline(x, y, fix_previous=False)
    p.add_func("monotone c2 (no fix)", f, interval, group='2m', color=p.color(6), **kwargs)
    p.add_func("monotone c2d1 (no fix)", f.derivative(), interval, 
               dash="dash", color=p.color(6,alpha=.5), group='2m', **kwargs)
    # Continuity with iterative monotone fix.
    f = monotone_quintic_spline(x,y)
    p.add_func("monotone c2", f, interval, group='2mf', color=p.color(7), **kwargs)
    p.add_func("monotone c2d1", f.derivative(), interval, 
               dash="dash", color=p.color(7,alpha=.5), group='2mf', **kwargs)

    p.show(file_name="monotone_quintic_interpolating_spline.html")
