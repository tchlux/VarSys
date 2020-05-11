# TODO: Local quintic fits can pass through nonmonotone pieces as long
#       as they meet the requirments. So make sure the first derivative
#       is zero at the points? Look into this, but it might cause
#       unwanted oscilation.


# This file provides the function `monotone_quintic_spline` which
# constructs spline interpolant over provided data. It starts by using
# divided differences to construct approximations for the first and
# second derivative at each of the data points. Then it iteratively
# adjusts the derivatives and second derivatives until all quintic
# polynomial pieces are monotone. A `Spline` object is returned.

from polynomial import Spline, fit, UsageError
IS_MONOTONE_CALLS = 0

# Convert a value or list of values into a list of Fraction objects.
def make_exact(value):
    from fraction import Fraction
    try:    return list(map(make_exact, value))
    except: return Fraction(value)

# Given x and y values, construct a monotonic quintic spline fit.
def monotone_quintic_spline(x, y=None, values=None, exact=True, 
                            max_steps=100, adaptive=1, verbose=False,
                            local_fits=None, free=False):
    from polynomial import Spline, local_facet, local_polynomial, polynomial
    if local_fits is None: local_fits = dict()
    if verbose:
        from polynomial import Polynomial
        print()
    # Convert 'x' to exact values.
    if exact: x = make_exact(x)

    # Build a set of 'values' if they were not provided.
    if (y is not None):
        if exact: y = make_exact(y)
        # Compute all known first and second derivative constraints.
        transitions = [i for i in range(1,len(y)-1)
                       if ((y[i]-y[i-1]) * (y[i+1]-y[i]) < 0)]
        flats = [i for i in range(1, len(y)-1)
                 if ((y[i]-y[i-1]) * (y[i+1]-y[i]) == 0)]
        # Build local quintic fits satisfying derivative constraints.
        values = [[y[i]] + [0]*2 for i in range(len(x))]
        for i in range(1,len(x)-1):
            left_slope = (y[i] - y[i-1]) / (x[i] - x[i-1])
            right_slope = (y[i+1] - y[i]) / (x[i+1] - x[i])
            if (left_slope * right_slope == 0): df = 0
            elif (left_slope * right_slope < 0) and (not free): df = 0
            else:
                idx = i
                if (i-1 in (flats + transitions)):
                    local_x, local_y = x[i-1:], y[i-1:]
                    idx = 1
                elif (i+1 in flats + transitions):
                    local_x, local_y = x[:i+2], y[:i+2]
                else:
                    local_x, local_y = x, y                    
                # # Local facet model.
                # f = local_facet(local_x, local_y, idx, size=3)

                # # Local quadratic interpolant.
                # f = local_polynomial(local_x, local_y, idx, order=3)

                # Local quadratic facet (minimum curvature).
                f1 = local_polynomial(local_x, local_y, idx-1, order=3)
                f2 = local_polynomial(local_x, local_y, idx,   order=3)
                f3 = local_polynomial(local_x, local_y, idx+1, order=3)
                ddf1 = abs(f1.derivative(2)(x[i]))
                ddf2 = abs(f2.derivative(2)(x[i]))
                ddf3 = abs(f3.derivative(2)(x[i]))
                if (ddf1 <= ddf2):
                    if (ddf1 <= ddf3): f = f1
                    else:              f = f3
                else:
                    if (ddf2 <= ddf3): f = f2
                    else:              f = f3

                # Compute the derivative.
                df = f.derivative()(x[i])

            # Store the derivative estimate.
            values[i][1] = df

        # Estimate the first derivative at the ends by using a local
        # quadratic that interpolates the first derivative at its neighbor.
        if (len(x) > 2):
            # Compute the left side.
            f = polynomial([x[0],x[1]], [y[0], y[1]], dx1=values[1][1])
            df = f.derivative()
            values[0][1:] = [df(x[0]), 0]
            local_fits[0] = (f, [x[0], x[1]])
            # Compute the right side.
            f = polynomial([x[-2],x[-1]], [y[-2], y[-1]], dx0=values[-2][1])
            df = f.derivative()
            values[-1][1:] = [df(x[-1]), 0]
            local_fits[len(x)-1] = (f, [x[-2], x[-1]])
        # Use a straight line inteprolant for two points.
        elif (len(x) == 2):
            values[0][1] = (y[1] - y[0]) / (x[1] - x[0])
            values[1][1] = values[0][1]
        print()

        # Estimate second derivatives.
        for i in range(1,len(x)-1):
            # Construct two quadratics that inteprolate this
            # derivative value and a neighbor, pick the one with the
            # minimum second derivative.
            f1 = polynomial([x[i-1],x[i]], [y[i-1],y[i]], dx1=values[i][1])
            f2 = polynomial([x[i],x[i+1]], [y[i],y[i+1]], dx0=values[i][1])
            ddf1 = f1.derivative(2)
            ddf2 = f2.derivative(2)
            if (abs(ddf1(x[i])) <= abs(ddf2(x[i]))): f = f1
            else:                                    f = f2
            # Compute the derivatives.
            ddf = f.derivative(2)
            # Store the values.
            values[i][2] = ddf(x[i])
            if verbose: print(f" estimate at {x[i]}: ", values[-1], Polynomial(f))

        print()

        # Store local fit information.
        print("Initial fits:")
        for i in range(len(x)):
            print("", i, "%.2f"%(x[i]),
                  "(%.2f %.2f %.2f)"%(y[i], values[i][1], values[i][2]))
            local_fits[i] = (polynomial([x[i]],[y[i]],
                                        dx0=values[i][1],
                                        ddx0=values[i][2]), 
                             [x[max(0,i-1)], x[min(len(x)-1,i+1)]])

        # Make all final values guesses exact.
        values = make_exact(values)

        # Print out the final set of values to the user.
        if verbose:
            print()
            print("Final values:")
            for i in range(len(x)):
                print(" at", x[i], "=", values[i])

    # If values were manually provided, use those.
    elif (values is not None):
        if exact: values = make_exact(values)
        # Verify monotonicity.
        if   all(values[i+1][0] >= values[i][0] for i in range(len(values)-1)): pass
        elif all(values[i][0] >= values[i+1][0] for i in range(len(values)-1)): pass
        else: raise(Exception("Provided 'y' data is not monotone."))
        # Make a spline over the points and derivative values.
        f = Spline(x, values)
        values = f.values
    else:
        raise(UsageError("Must provided either a flat 'y' array or a (N,3) values array."))

    # Make all pieces monotone.
    initial_values = [v.copy() for v in values]
    check_counts = {}
    change_counts = {}
    change_values = {i:values[i][1] / max_steps for i in range(len(values))}
    to_check = list(range(len(values)-1))
    if verbose:
        print()
        print("Starting quintic monotonicity fix, smallest allowed steps:")
        for k in sorted(change_values):
            print(f" at {k} =",change_values[k])
        print()
    # Cycle over all intervals that need to be checked for monotonicity.
    while (len(to_check) > 0):
        i = to_check.pop(0)
        if (values[i][1] * values[i+1][1] < 0) and free: continue
        if verbose: print(f"\n  checking ({i+1}) ..")
        # Update this interval to make it monotone.
        initial_values = [values[i].copy(), values[i+1].copy()]
        changed = monotone_quintic(
            [x[i], x[i+1]], [values[i], values[i+1]])
        check_counts[i] = check_counts.get(i,0) + 1
        # If in fact this interval wasn't changed, remove the increment.
        if (not changed) and verbose: print(f"   not changed..")
        # Mark the adjacent intervals as needing to be checked,
        # since this adjustment may have broken their monotonicity.
        if changed:
            # Track the number of times this particular interval has been changed.
            change_counts[i] = change_counts.get(i,0) + 1
            # Shrink neighboring derivative values if this interval is looping.
            if change_counts[i] > 1:
                # Compute the amount by which the second derivative
                # had violated the monotonicity constraint.
                for j, init_vals in zip([i,i+1], initial_values):
                    # Skip derivatives that have already been zeroed.
                    if values[j][1] == 0: continue
                    # Compute the second derivative change (violation)
                    # and the magnitude of the derivative at each end.
                    violation = abs(init_vals[2] - values[j][2])
                    scale = abs(values[j][1])
                    # Compute the relative step size using a linear
                    # fractional transformation.
                    relative_change = 1 - 1 / (1 + adaptive*(violation + 1/scale))
                    change = max(change_values[j], values[j][1] * relative_change)
                    # Make sure the change does not push past 0.
                    change = min(change, values[j][1])
                    values[j][1] -= change
                change_counts[i] = 0
                change_counts[i+1] = 0
            # Queue up this interval and its neighbors to be checked again.
            next_value  = to_check[0]  if (len(to_check) > 0) else None
            second_last = to_check[-2] if (len(to_check) > 1) else None
            last_value  = to_check[-1] if (len(to_check) > 0) else None
            if (i > 0) and (i-1 not in {second_last, last_value}):
                to_check.append( i-1 )
            if (i not in {last_value}):
                to_check.append( i )
            if (i+1 < len(to_check)) and (i+1 != next_value):
                to_check.insert(0, i+1)
            # Maintain exactness if necessary!
            if exact:
                values[i][:] = make_exact(values[i])
                values[i+1][:] = make_exact(values[i+1])
            # Show some printouts to the user for verbosity.
            if verbose:
                print(f"   [update {change_counts[i]}] interval {i+1} changed ({changed})..")
                print(f"     {float(x[i]):.3f} {float(x[i+1]):.3f} ({float(values[i][0]):.2e} {float(values[i+1][0]):.2e})")
                print( "     ", ' '.join([f"{float(v): .5e}" for v in values[i]]))
                print( "     ", ' '.join([f"{float(v): .5e}" for v in values[i+1]]))
        
    print()
    print("Final fits:")
    for i in range(len(x)):
        print("", i, "%.2f"%(x[i]),
              "(%.2f %.2f %.2f)"%(y[i], values[i][1], values[i][2]))

    # Create "smooth estimates of all points"

    # # Manually return the value to its original.
    # values[4][2] = -.67

    # Construct a new spline over the (updated) values for derivatives.
    return Spline(x, values)


# Generate a minimum curvature estimate of Df and DDf at x[i] by using
# neighboring values and second derivatives.
def minimum_curvature_estimate(x, fx, i, use=1):
    if (len(x) == 3):
    f1 = polynomial([x[i-1],x[i]], [y[i-1],y[i]], dx1=values[i][1])
    f2 = polynomial([x[i],x[i+1]], [y[i],y[i+1]], dx0=values[i][1])
    ddf1 = f1.derivative(2)
    ddf2 = f2.derivative(2)
    if (abs(ddf1(x[i])) <= abs(ddf2(x[i]))): f = f1
    else:                                    f = f2
    # Compute the derivatives.
    ddf = f.derivative(2)

    


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
        changed = any(v != 0 for v in (DX0, DX1, DDX0, DDX1)) * 1
        y[0][:] = X0, 0, 0
        y[1][:] = X1, 0, 0
        return changed
    # Convert all the tests to the monotone increasing case (will undo).
    sign = (-1) ** int(function_change < 0)
    X0 *= sign
    X1 *= sign
    DX0 *= sign
    DX1 *= sign
    DDX0 *= sign
    DDX1 *= sign
    interval_slope *= sign
    # Set DX0 and DX1 to be the median of these three choices.
    if not (0 <= DX0 <= 14*interval_slope):
        new_val = (sorted([0, DX0, 14*interval_slope])[1])
        changed = (new_val != DX0)
        if changed: DX0 = new_val
    if not (0 <= DX1 <= 14*interval_slope):
        new_val = (sorted([0, DX1, 14*interval_slope])[1])
        changed += (new_val != DX1) + 10**1
        if changed: DX1 = new_val
    # Compute repeatedly used values "A" (left ratio) and "B" (right ratio).
    A = DX0 / interval_slope
    B = DX1 / interval_slope
    assert(A >= 0)
    assert(B >= 0)
    # Use a (simplified) monotone cubic over this region if AB = 0.
    # Only consider the coefficients less than the x^4, because that
    # term is strictly positive.
    if (A*B < accuracy):
        # print("interval: (%.2f   %.2f)"%(U0,U1))
        # Compute a temporary variable that is useful in the algebra.
        w = U0 - U1
        # First, find a DX0 and DX1 that makes DDX0 have nonempty feasible region.
        DX_multiplier = max(0, (- 20*(X1-X0) / w) / (5*DX0 + 4*DX1))
        if (DX_multiplier < 1):
            DX0 *= DX_multiplier
            DX1 *= DX_multiplier
        # Second, cap DDX0 so that the DDX1 feasible region is nonempty.
        max_DDX0 = (4*(2*DX0 + DX1) - 20*(X0-X1)/w) / w
        # print("   max_DDX0: %.2f"%(max_DDX0))
        if (DDX0 > max_DDX0):
            DDX0 = max_DDX0
        # Enforce \gamma >= \delta
        min_DDX0 = 3 * DX0 / w
        # print("   min_DDX0: %.2f"%(min_DDX0))
        if (min_DDX0 > DDX0):
            DDX0 = min_DDX0
        # Enforce \alpha >= 0 
        max_DDX1 = -4*DX1 / w
        # print("   max_DDX1: %.2f"%(max_DDX1))
        if (DDX1 > max_DDX1):
            DDX1 = max_DDX1
        # Enforce \beta >= \alpha
        min_DDX1 = (3*DDX0*w - (24*DX0 + 32*DX1) + 60*(X0-X1)/w) / (5*w)
        # print("   min_DDX1: %.2f"%(min_DDX1))
        if (DDX1 < min_DDX1):
            DDX1 = min_DDX1
        # Check for contradictions, which should never happen!
        assert(min_DDX0 - max_DDX0 <= 2**(-26))
        assert(min_DDX1 - max_DDX1 <= 2**(-26))

        # Function for telling if the given values are monotone.
        def is_monotone():
            alpha = (4*DX1 + DDX1*(U0-U1)) * (U0-U1) / (X0-X1)
            beta = 30 + ((-24*(DX0+DX1) + 3*(DDX0-DDX1)*(U0-U1))*(U0-U1)) / (2 * (X0-X1))
            gamma = ((U0-U1) * (4*DX0 + DDX0*(U1-U0))) / (X0-X1)
            delta = DX0*(U0-U1) / (X0-X1)
            return (alpha >= 0) and (delta >= 0) and \
                (beta >= alpha - (4*alpha*delta)**(1/2)) and \
                (gamma >= delta - (4*alpha*delta)**(1/2))

        assert(is_monotone())
        # Known monotone values for DDX0 and DDX1
        target_DX0 = DX0
        target_DX1 = DX1
        target_DDX0 = DDX0
        target_DDX1 = DDX1
        # Reset to original values (sign adjusted).
        DX0, DDX0 = y[0][1:]
        DX1, DDX1 = y[1][1:]
        DX0 *= sign
        DX1 *= sign
        DDX0 *= sign
        DDX1 *= sign
        # Do a binary search between known feasible values and current values.
        if not (is_monotone()):
            original_DX0, original_DDX0 = DX0, DDX0
            original_DX1, original_DDX1 = DX1, DDX1
            # If this function is not monotone, perform a binary
            # search for the nearest-to-original monotone DDX values.
            low_bound, upp_bound = 0, 1
            # Continue dividing the interval in 2 to find the smallest
            # "upper bound" (amount target) that satisfies monotonicity.
            while ((upp_bound - low_bound) > accuracy):
                to_target = upp_bound / 2  +  low_bound / 2
                # If we found the limit of floating point numbers, break.
                if ((to_target == upp_bound) or (to_target == low_bound)): break
                # Recompute DDX0 and DDX1 based on the to_target.
                DX0 = (1-to_target) * original_DX0 + to_target * target_DX0
                DX1 = (1-to_target) * original_DX1 + to_target * target_DX1
                DDX0 = (1-to_target) * original_DDX0 + to_target * target_DDX0
                DDX1 = (1-to_target) * original_DDX1 + to_target * target_DDX1
                # Otherwise, proceed with a binary seaarch.
                if is_monotone(): upp_bound = to_target
                else:             low_bound = to_target
            # Store the smallest amount "target" for DDX0 and DDX1 that is monotone.
            DX0 = (1-upp_bound) * original_DX0 + upp_bound * target_DX0
            DX1 = (1-upp_bound) * original_DX1 + upp_bound * target_DX1
            DDX0 = (1-upp_bound) * original_DDX0 + upp_bound * target_DDX0
            DDX1 = (1-upp_bound) * original_DDX1 + upp_bound * target_DDX1
            assert(is_monotone())

        # If any values have changed, record that.
        if not all(old==new for (old,new) in zip(
                y[0][1:]+y[1][1:], [sign*DX0,sign*DDX0,sign*DX1,sign*DDX1])):
            changed = True
            y[0][1:] = sign*DX0, sign*DDX0
            y[1][1:] = sign*DX1, sign*DDX1

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
        global IS_MONOTONE_CALLS
        IS_MONOTONE_CALLS += 1
        a = alpha_constant + alpha_multiplier * DDX1
        g = gamma_constant + gamma_multiplier * DDX0
        b = beta_constant  + beta_multiplier  * (DDX0 - DDX1)
        if b <= 6: bound = - (b + 2) / 2
        else:      bound = -2 * (b - 2)**(1/2)
        return (a-accuracy > bound) and (g-accuracy > bound)

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
    y[0][:] = sign*X0, sign*DX0, sign*DDX0
    y[1][:] = sign*X1, sign*DX1, sign*DDX1
    return changed


# Given x and y values, construct a monotonic quintic spline fit.
def monotone_cubic_spline(x, y=None, values=None, exact=False, **kwargs):
    import numpy as np
    from polynomial import Spline
    # Convert 'x' to exact values.
    if exact: x = make_exact(x)
    # Process the provided data to prepare for monotonicity checks.
    if (y is not None):
        if exact: y = make_exact(y)
        
        kwarg = kwargs.copy()
        if   all(y[i+1] >= y[i] for i in range(len(y)-1)): kwarg.update(dict(min_d1=0))
        elif all(y[i] >= y[i+1] for i in range(len(y)-1)): kwarg.update(dict(max_d1=0))
        else: raise(Exception("Provided 'y' data is not monotone."))
        f = fit(x, y, continuity=1, **kwarg)
        values = f.values
        # # Make the first derivative 0 at the right end.
        # values[-1][1] = 0
    elif (values is not None):
        if exact: values = make_exact(values)
        # Verify monotonicity.
        if   all(values[i+1][0] >= values[i][0] for i in range(len(values)-1)): pass
        elif all(values[i][0] >= values[i+1][0] for i in range(len(values)-1)): pass
        else: raise(Exception("Provided 'y' data is not monotone."))
        # Make a spline over the points and derivative values.
        f = Spline(x, values)
        values = f.values
    else:
        raise(UsageError("Must provided either a flat 'y' array or a (N,2) values array."))
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
    # 1, 4, (None)      -> Bad (now good) (not monotone) [bad first derivative conditions]
    # 0, 13, (None)     -> Bad (now good) (not monotone after first pass) [no previous-interval fix]
    # 0, 30, (-3,None)  -> Bad (now good) (far right still not monotone after passes) [bad simplied conditions]
    # 0, 100, (None)    -> Ehh (now good) (barely non-monotone after 100 steps) [made lower ceiling for derivatives]
    # 0, 1000, (None)   -> Good
    # 
    # --------------------------------------------------------------------

    # SEED = 0
    # NODES = 100
    # SUBSET = slice(None) 

    SEED = 9
    INTERNAL_NODES = 10
    SUBSET = slice(None)

    # Generate random data to test the monotonic fit function.
    import numpy as np
    np.random.seed(SEED)
    nodes = INTERNAL_NODES + 2
    x = np.linspace(0, 1, nodes)
    y = sorted(np.random.normal(size=(nodes,)))
    x -= min(x); x /= max(x)
    y -= min(y); y /= max(y); y -= max(y)
    # Convert these arrays to lists.
    x, y = list(x)[SUBSET], list(y)[SUBSET]

    # Convert these arrays to exact arithmetic.
    x = make_exact(x)
    y = make_exact(y)
    interval = [float(min(x)), float(max(x))]

    # Generate a plot to see what it all looks like.
    from util.plot import Plot
    # Create a holder for the plot.
    p = Plot()
    p.add("Points", list(map(float, x)), list(map(float, y)))
    kwargs = dict(plot_points=1000)
    # Continuity 0
    p.add_func("continuity 0", fit(x,y,continuity=0),
               interval, group='0', **kwargs)
    p.add_func("c0 d1", fit(x,y,continuity=0).derivative(),
               interval, dash="dash", color=p.color(p.color_num,alpha=.5), 
               group='0', **kwargs)
    # Continuity 1
    f = monotone_cubic_spline(x,y)
    p.add_func("monotone c1", f, interval, group='1m', **kwargs)
    p.add_func("monotone c1 d1", f.derivative(), interval, 
               dash="dash", color=p.color(p.color_num,alpha=.5), group='1m', **kwargs)
    # Continuity 2
    f,_,_,_ = monotone_quintic_spline(x,y, verbose=True)
    # f = monotone_quintic_spline(x,y)
    kwargs = dict(plot_points=10000)
    p.add_func("monotone c2", f, interval, group='2mf', color=p.color(7), **kwargs)
    p.add_func("monotone c2d1", f.derivative(), interval, 
               dash="dash", color=p.color(7,alpha=.5), group='2mf', **kwargs)
    # Show the plot.
    p.show(file_name="monotone_quintic_interpolating_spline.html")
