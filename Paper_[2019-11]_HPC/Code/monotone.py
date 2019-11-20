# This file provides the function `monotone_quintic_spline` which
# constructs spline interpolant over provided data. It starts by using
# divided differences to construct approximations for the first and
# second derivative at each of the data points. Then it iteratively
# adjusts the derivatives and second derivatives until all quintic
# polynomial pieces are monotone. A `Spline` object is returned.

from polynomial import fit

# Given x and y values, construct a monotonic quintic spline fit.
def monotone_quintic_spline(x, y, ends=2, mids=2, exact=True,
                            max_steps=100, verbose=False):
    if verbose: print()
    from polynomial import Spline
    if exact:
        from util.math import Fraction
        x = list(map(Fraction, x))
        y = list(map(Fraction, y))
    # Check for which type of monotonicity this function maintains.
    kwarg = {}
    if   all(y[i+1] >= y[i] for i in range(len(y)-1)): kwarg = dict(min_d1=0)
    elif all(y[i] >= y[i+1] for i in range(len(y)-1)): kwarg = dict(max_d1=0)
    else: raise(Exception("Provided 'y' data is not monotone."))
    # Construct an initial fit that is twice continuous.
    f = fit(x, y, continuity=2, ends=ends, mids=mids, **kwarg)
    values = f.values
    # Make all pieces monotone.
    change_counts = {}
    to_check = list(range(len(values)-1))
    if verbose: print("Starting quintic monotonicity fix..\n")
    # Cycle over all intervals that need to be checked for monotonicity.
    while (len(to_check) > 0):
        i = to_check.pop(0)
        # Track the number of times this particular interval has been checked.
        change_counts[i] = change_counts.get(i,0) + 1
        # Do NOT update this interval if it has been changed too many times.
        step = change_counts[i]
        if (step > max_steps): continue
        if verbose: print(f"  checking ({i+1}) ..")
        # Update this interval to make it monotone.
        changed = monotone_quintic(
            [x[i], x[i+1]], [values[i], values[i+1]])
        # Mark the adjacent intervals as needing to be checked,
        # since this adjustment may have broken their monotonicity.
        if changed:
            # Queue up this interval and its neighbors to be checked again.
            if (i > 0):               to_check.append( i-1 )
            if (step < max_steps):    to_check.append(i)
            if (i+1 < len(to_check)): to_check.insert(0, i+1)
            # Maintain exactness if necessary!
            if exact:
                values[i][:] = map(Fraction, values[i])
                values[i+1][:] = map(Fraction, values[i+1])
            # Show some printouts to the user for verbosity.
            if verbose:
                print(f"   [update {step}] interval {i+1} changed..")
                print(f"     {float(x[i]):.2f} {float(x[i+1]):.2f}")
                print( "     ", ' '.join([f"{float(v): .2f}" for v in values[i]]))
                print( "     ", ' '.join([f"{float(v): .2f}" for v in values[i+1]]))
    if verbose: print("done.\n")
    # Construct a new spline over the (updated) values for derivatives.
    return Spline(f.knots, f.values)

# Given a (x1, x2) and ([y1, d1y1, d2y1], [y2, d1y2, d2y2]), compute
# the rescaled derivative values for a monotone quintic piece.
def monotone_quintic(x, y, accuracy=2**(-26)):
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
    # Use a (simplified) monotone cubic over this region if AB = 0.
    # Only consider the coefficients less than the x^4, because that
    # term is strictly positive.
    if (A*B == 0):
        # Convert all the tests to the monotone increasing case (will undo).
        X0 *= sign
        X1 *= sign
        DX0 *= sign
        DX1 *= sign
        DDX0 *= sign
        DDX1 *= sign
        # Construction functions for checking the feasible DDX1 region.
        w = U0 - U1
        max_DDX1 = lambda: -4*DX1 / w
        min_DDX1 = lambda: (60*(X0-X1)/w - (24*DX0 + 32*DX1)) / (5*w)
        # Store original values for reverting after searches.
        orig_DX1 = DX1
        orig_DX0 = DX0
        # Make sure that given a viable DX0 and DDX0, there is an achievable DX1 (>= 0).
        DX1 = 0
        #   if there is currently not an achievable solution, change DX0.
        if (max_DDX1() < min_DDX1()):
            low = 0
            upp = DX0
            while ((upp - low) > accuracy):
                DX0 = (low+upp) / 2
                if (max_DDX1() < min_DDX1()): upp = DX0
                else:                         low = DX0
            DX0 = low
        #   update the new first derivative on left if it changed.
        if (orig_DX0 != DX0):
            y[0][1] = DX0*sign
            changed = True
        # Make sure that given DX1, there is an achievable DDX1.
        DX1 = orig_DX1
        #   if there is currently not an achievable solutoin, change DX1.
        if (max_DDX1() < min_DDX1()):
            low = 0
            upp = DX1
            while ((upp - low) > accuracy):
                DX1 = (low+upp) / 2
                if (max_DDX1() < min_DDX1()): upp = DX1
                else:                         low = DX1
            DX1 = low
        #   update the new first derivative on right if it changed.
        if (orig_DX1 != DX1):
            y[1][1] = DX1*sign
            changed = True
        # Enforce \alpha >= 0 
        max_DDX1 = -4*DX1 / w
        if (DDX1 > max_DDX1):
            DDX1 = max_DDX1
            y[1][2] = DDX1*sign
            changed = True
        # Enforce \beta >= \alpha
        min_DDX1 = (3*DDX0*w - (24*DX0 + 32*DX1) + (60*(X0-X1))/w) / (5*w)
        if (DDX1 < min_DDX1):
            DDX1 = min_DDX1
            y[1][2] = DDX1*sign
            changed = True
        # Enforce \gamma >= \delta
        min_DDX0 = 3 * DX0 / w
        if (min_DDX0 > DDX0):
            DDX0 = min_DDX0
            y[0][2] = DDX0*sign
            changed = True
        # Return whether or not the values were changed.
        return changed

    # Clip derivative values that are too large (to ensure that
    # shrinking the derivative vectors on either end will not break
    # monotonicity). (clipping at 6 box is enough, with 8 needs more treatment)
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

    # Update "y" and return the updated version.
    y[0][:] = X0, DX0, DDX0
    y[1][:] = X1, DX1, DDX1
    return changed


# Given x and y values, construct a monotonic quintic spline fit.
def monotone_cubic_spline(x, y):
    from polynomial import Spline
    kwarg = {}
    if   all(y[i+1] >= y[i] for i in range(len(y)-1)): kwarg = dict(min_d1=0)
    elif all(y[i] >= y[i+1] for i in range(len(y)-1)): kwarg = dict(max_d1=0)
    else: raise(Exception("Provided 'y' data is not monotone."))
    f = fit(x, y, continuity=1, **kwarg)
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



if __name__ == "__main__":
    # --------------------------------------------------------------------
    #               TEST CASES
    # 
    # 0, 4, (None)      -> Good
    # 1, 4, (None)      -> Bad (not monotone) (now good) [bad first derivative conditions]
    # 0, 13, (None)     -> Bad (not monotone after first pass) (now good) [no previous-interval fix]
    # 0, 30, (-3,None)  -> Bad (far right still not monotone after passes) (now good) [bad simplied conditions]
    # 0, 100, (None)    -> Ehh (some intervals are barely non-monotone after 100 steps)
    # 0, 1000, (None)   -> Ehh (same as above)
    # 
    # --------------------------------------------------------------------

    # SEED = 0
    # NODES = 100
    # SUBSET = slice(None) 

    SEED = 0
    NODES = 13
    SUBSET = slice(None)

    # Generate random data to test the monotonic fit function.
    import numpy as np
    np.random.seed(SEED)
    nodes = NODES + 2
    x = np.linspace(0, 1, nodes)
    y = sorted(np.random.normal(size=(nodes,)))
    x -= min(x); x /= max(x)
    y -= min(y); y /= max(y); y -= max(y)
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
    p.add_func("continuity 0", fit(x,y,continuity=0),
               interval, group='0', **kwargs)
    p.add_func("c0 d1", fit(x,y,continuity=0, min_d1=0).derivative(),
               interval, dash="dash", color=p.color(p.color_num,alpha=.5), 
               group='0', **kwargs)
    # Continuity 1
    # p.add_func("continuity 1", fit(x,y,continuity=1, ends=1, mids=0), 
    #            interval, group='1', **kwargs)
    # p.add_func("c1 d1", fit(x,y,continuity=1, min_d1=0).derivative(), 
    #            interval, dash="dash", color=p.color(p.color_num,alpha=.5), 
    #            group='1', **kwargs)
    f = monotone_cubic_spline(x,y)
    p.add_func("monotone c1", f, interval, group='1m', **kwargs)
    p.add_func("monotone c1 d1", f.derivative(), interval, 
               dash="dash", color=p.color(p.color_num,alpha=.5), group='1m', **kwargs)
    # Continuity 2
    # p.add_func("continuity 2", fit(x,y,continuity=2, min_d1=0), 
    #            interval, group='2', **kwargs)
    # p.add_func("c2 d1", fit(x,y,continuity=2, min_d1=0).derivative(), 
    #            interval, dash="dash", color=p.color(p.color_num,alpha=.5), group='2', **kwargs)
    print("Constructing monotone quintic..")
    f = monotone_quintic_spline(x,y, verbose=True)
    print("done.")
    # kwargs = dict(plot_points=10000)
    p.add_func("monotone c2", f, interval, group='2mf', color=p.color(7), **kwargs)
    p.add_func("monotone c2d1", f.derivative(), interval, 
               dash="dash", color=p.color(7,alpha=.5), group='2mf', **kwargs)

    p.show(file_name="monotone_quintic_interpolating_spline.html")
