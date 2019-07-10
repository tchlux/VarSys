# - implement the four solutions that can replace the binary search
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
def monotone_quintic_spline(x, y, ends=2, mids=2, fix_previous=True):
    from polynomial import Spline
    f = fit(x, y, continuity=2, non_decreasing=True, ends=ends, mids=mids)
    values = f.values
    # Make all pieces monotone.
    i = 0
    while i < len(values)-1:
        changed = monotone_quintic(
            [x[i], x[i+1]], [values[i], values[i+1]])
        if fix_previous and changed and (i > 0):
            step = 0
            changed = False
            # While the previous interval is also broken,
            while not changed:
                step += 1
                # Shrink the derivative value at this point.
                values[i][1] *= .99

                # Fix the previous interval.
                changed = monotone_quintic(
                    [x[i-1], x[i]], [values[i-1], values[i]])
                # Fix this interval again.
                changed = monotone_quintic(
                    [x[i], x[i+1]], [values[i], values[i+1]]) or changed
                if step >= 1000:
                    print()
                    print('-'*70)
                    print()
                    print(f"Interval {i}: {x[i]:.3f} {x[i+1]:.3f}")
                    print()
                    print('-'*70)
                    print()
                    raise(Exception("Error for this interval.."))
            print("Steps:", step)
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
        monotone_cubic([x[i], x[i+1]], [values[i], values[i+1]])
    # Construct a new spline over the (updated) values for derivatives.
    return Spline(f.knots, f.values)

# Given a (x1, x2) and ([y1, d1y1], [y2, d1y2]), compute the rescaled
# y values for a monotone cubic piece.
def monotone_cubic(x,y):
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
    sign = (-1) ** int(function_change < 0)
    # Set DX0 and DX1 to be the median of these three choices.
    if not (0 <= sign*DX0 <= sign*14*interval_slope):
        DX0 = sign * (sorted([0, sign*DX0, sign*14*interval_slope])[1])
        changed = True
    if not (0 <= sign*DX1 <= sign*14*interval_slope):
        DX1 = sign * (sorted([0, sign*DX1, sign*14*interval_slope])[1])
        changed = True
    # Compute repeatedly used values "A" (left ratio) and "B" (right ratio).
    A = DX0 / interval_slope
    B = DX1 / interval_slope
    assert(A >= 0)
    assert(B >= 0)
    # Use a monotone cubic over this region if AB = 0.
    if (A*B <= 0):
        y[0][:] = X0, DX0, DDX0
        y[1][:] = X1, DX1, DDX1
        return monotone_cubic(x, y) or True

    # # Scale the derivative vector to make tau_1 positive.
    # tau_1 = 24 + 2*(A*B)**(1/2) - 3*(A+B)
    # if (tau_1 <= 0):
    #     # Compute the rescale factor necessary to make tau_1 0.
    #     rescale = 24 * (3*(A+B) + 2*(A*B)**(1/2)) / (9*(A**2+B**2) + 14*(A*B))
    #     rescale -= 2**(-26) * rescale
    #     # Rescale the derivatives
    #     DX0 *= rescale
    #     DX1 *= rescale
    #     # Recompute A and B.
    #     A = DX0 / interval_slope
    #     B = DX1 / interval_slope
    #     # Record the change.
    #     changed = True

    # Clip derivative values that are too large (to ensure that
    # shrinking the derivative vectors on either end will not break
    # monotonicity). (clipping at 6 alone is enough, with 8 needs more)
    mult = 6 / max(A, B)
    if (mult < 1):
        DX0 *= mult
        DX1 *= mult
        A = DX0 / interval_slope
        B = DX1 / interval_slope
        changed = True
    # Make sure that the first monotonicity condition is met.
    tau_1 = 24 + 2*(A*B)**(1/2) - 3*(A+B)
    try: assert(tau_1 > 0)
    except:
        class BadTau(Exception): pass
        raise(NonMonotoneTau(f"Bad Tau 1 value: {tau_1}"))

    # Compute DDX0 and DDX1 that satisfy monotonicity by scaling (C,D)
    # down until montonicity is achieved (using binary search).
    alpha_constant   = 4 * (B**(1/4) / A**(1/4))
    alpha_multiplier = ((U0-U1) / DX1) * B**(1/4) / A**(1/4)
    gamma_constant   = 4 * ((DX0 * B**(3/4)) / (DX1 * A**(3/4)))
    gamma_multiplier = ((U1-U0) / DX1) * B**(3/4) / A**(3/4)
    beta_constant    = (12 * (DX0+DX1) * (U1-U0) + 30 * (X0-X1)) / ((X0-X1) * A**(1/2) * B**(1/2))
    beta_multiplier  = (3 * (U0-U1)**2) / (2 * (X0-X1) * A**(1/2) * B**(1/2))
    def is_monotone():
        a = alpha_constant + alpha_multiplier * DDX1
        g = gamma_constant + gamma_multiplier * DDX0
        b = beta_constant  + beta_multiplier  * (DDX0 - DDX1)
        if b <= 6: bound = - (b + 2) / 2
        else:      bound = -2 * (b - 2)**(1/2)
        return (a > bound) and (g > bound)
    def print_abg():
        DDX0 = ratio * original_DDX0 + (1-ratio) * target_DDX0
        DDX1 = ratio * original_DDX1 + (1-ratio) * target_DDX1
        a = alpha_constant + alpha_multiplier * DDX1
        g = gamma_constant + gamma_multiplier * DDX0
        b = beta_constant  + beta_multiplier  * (DDX0 - DDX1)
        print()
        print(f"ratio = {ratio}")
        print(f"DDX0 = {DDX0}")
        print(f"DDX1 = {DDX1}")
        print(f" alpha = {a}")
        print(f" gamma = {g}")
        print(f" beta  = {b}")
        if b <= 6:
            bound = - (b + 2) / 2
            print("  SMALL, b <= 6")
        else:
            bound = -2 * (b - 2)**(1/2)
            print("  BIG,   b > 6")
        print(f"  bound = {bound}")
        print(f"  alpha > bound  ->  {a > bound}")
        print(f"  gamma > bound  ->  {g > bound}")
        print()

    # Move the second derivative towards a working value until
    # monotonicity is achieved.
    target_DDX0 = - A**(1/2) * (7*A**(1/2) + 3*B**(1/2)) * interval_slope / interval_width
    target_DDX1 =   B**(1/2) * (3*A**(1/2) + 7*B**(1/2)) * interval_slope / interval_width

    # If this function is not monotone, perform a binary
    # search for the nearest-to-original monotone DDX values.
    if not is_monotone():
        accuracy = 2**(-26) # = SQRT(EPSILON( 0.0_REAL64 ))
        original_DDX0, original_DDX1 = DDX0, DDX1

        print("-" * 70)
        print()
        print("DDX0:  ",float(original_DDX0))
        print("DDX1:  ",float(original_DDX1))
        print()
        print("DDX TARGET")

        DDX0 = target_DDX0
        DDX1 = target_DDX1
        ratio = 0
        print_abg()

        from math import sqrt
        # CASE 1:  (alpha == bound) and (beta <= 6)
        try:
            print("CASE 1   (alpha == bound) and  (beta <= 6)")
            ratio = (3*B**0.25*DX0*(U0 - U1)**2*(7*DX0*(U0 - U1) + 7*DX1*(U0 - U1) + 4*(original_DDX0 - original_DDX1)*(X0 - X1)) + 2*sqrt(A)*B**0.75*DX0*(9*U0**2 - 18*U0*U1 + 9*U1**2 + 8*X0 - 8*X1)*(X0 - X1) - 12*A**1.75*sqrt(B)*(U0 - U1)*(X0 - X1)**2 - 4*A**1.25*(U0 - U1)*(X0 - X1)*(7*DX1*(U0 - U1) + 4*original_DDX1*(-X0 + X1)))/(3*B**0.25*DX0*(DX0*(U0 - U1)*(7*U0**2 - 14*U0*U1 + 7*U1**2 + 32*X0 - 32*X1) + DX1*(U0 - U1)*(7*U0**2 - 14*U0*U1 + 7*U1**2 + 32*X0 - 32*X1) - 80*(X0 - X1)**2) + 18*sqrt(A)*B**0.75*DX0*(U0 - U1)**2*(X0 - X1) - 4*A**1.25*DX1*(7*U0**2 - 14*U0*U1 + 7*U1**2 + 16*(X0 - X1))*(X0 - X1) - 12*A**1.75*sqrt(B)*(U0 - U1)*(X0 - X1)**2)
            print_abg()
        except:
            print()
            print("Failed..")
            print()
        # CASE 2:  (gamma == bound) and (beta <= 6)
        try:
            print("CASE 2   (gamma == bound) and  (beta <= 6)")
            ratio = (21*DX0**2*(U0 - U1)**3 + 21*DX0*DX1*(U0 - U1)**3 + 2*DX0*(-14*A**0.75*B**0.25*(U0 - U1)**2 + 6*(original_DDX0 - original_DDX1)*(U0 - U1)**2 + sqrt(A)*sqrt(B)*(9*U0**2 - 18*U0*U1 + 9*U1**2 + 8*X0 - 8*X1))*(X0 - X1) - 4*A**0.75*B**0.25*(3*sqrt(A)*sqrt(B) + 4*original_DDX0)*(U0 - U1)*(X0 - X1)**2)/(3*DX0**2*(U0 - U1)*(7*U0**2 - 14*U0*U1 + 7*U1**2 + 32*(X0 - X1)) + 3*DX0*DX1*(U0 - U1)*(7*U0**2 - 14*U0*U1 + 7*U1**2 + 32*(X0 - X1)) - 2*DX0*(-9*sqrt(A)*sqrt(B)*(U0 - U1)**2 + 2*A**0.75*B**0.25*(7*U0**2 - 14*U0*U1 + 7*U1**2 + 16*X0 - 16*X1) + 120*(X0 - X1))*(X0 - X1) - 12*A**1.25*B**0.75*(U0 - U1)*(X0 - X1)**2)
            print_abg()
        except:
            print()
            print("Failed..")
            print()
        # CASE 3:  (alpha == bound) and (beta > 6)
        try:
            print("CASE 3a  (alpha == bound) and (beta > 6)")
            ratio = (B**1.5*DX0**2*(X0 - X1)**2*((A**1.5*(U0 - U1)*(7*DX1*(U0 - U1) + (3*sqrt(A)*sqrt(B) - 4*original_DDX1)*(X0 - X1))*(DX1*(7*U0**2 - 14*U0*U1 + 7*U1**2 + 16*X0 - 16*X1) + 3*sqrt(A)*sqrt(B)*(U0 - U1)*(X0 - X1)))/(B**1.5*DX0**2*(X0 - X1)**2) - (12*(DX0*(U0 - U1)*(7*U0**2 - 14*U0*U1 + 7*U1**2 + 32*(X0 - X1)) + DX1*(U0 - U1)*(7*U0**2 - 14*U0*U1 + 7*U1**2 + 32*(X0 - X1)) + 2*(X0 - X1)*(3*sqrt(A)*sqrt(B)*(U0 - U1)**2 + 40*(-X0 + X1))))/(sqrt(A)*sqrt(B)*(X0 - X1)**2) - 2*sqrt(2)*sqrt((3*(U0 - U1)**2*(7*DX0*(U0 - U1) + 7*DX1*(U0 - U1) + 2*(3*sqrt(A)*sqrt(B) + 2*original_DDX0 - 2*original_DDX1)*(X0 - X1))*(A*DX1*(7*U0**2 - 14*U0*U1 + 7*U1**2 + 16*X0 - 16*X1) + 3*A**1.5*sqrt(B)*(U0 - U1)*(X0 - X1))**2 - 16*A**2.5*sqrt(B)*(DX1*(7*U0**2 - 14*U0*U1 + 7*U1**2 + 16*X0 - 16*X1) + 3*sqrt(A)*sqrt(B)*(U0 - U1)*(X0 - X1))**2*(X0 - X1)**2 - 3*A**2*(U0 - U1)*(7*DX1*(U0 - U1) + (3*sqrt(A)*sqrt(B) - 4*original_DDX1)*(X0 - X1))*(DX1*(7*U0**2 - 14*U0*U1 + 7*U1**2 + 16*X0 - 16*X1) + 3*sqrt(A)*sqrt(B)*(U0 - U1)*(X0 - X1))*(DX0*(U0 - U1)*(7*U0**2 - 14*U0*U1 + 7*U1**2 + 32*(X0 - X1)) + DX1*(U0 - U1)*(7*U0**2 - 14*U0*U1 + 7*U1**2 + 32*(X0 - X1)) + 2*(X0 - X1)*(3*sqrt(A)*sqrt(B)*(U0 - U1)**2 + 40*(-X0 + X1))) + 18*B*DX0**2*(DX0*(U0 - U1)*(7*U0**2 - 14*U0*U1 + 7*U1**2 + 32*(X0 - X1)) + DX1*(U0 - U1)*(7*U0**2 - 14*U0*U1 + 7*U1**2 + 32*(X0 - X1)) + 2*(X0 - X1)*(3*sqrt(A)*sqrt(B)*(U0 - U1)**2 + 40*(-X0 + X1)))**2)/(A*B**2*DX0**2*(X0 - X1)**4))))/(A**1.5*(DX1*(7*U0**2 - 14*U0*U1 + 7*U1**2 + 16*X0 - 16*X1) + 3*sqrt(A)*sqrt(B)*(U0 - U1)*(X0 - X1))**2)
            print_abg()
        except:
            print()
            print("Failed..")
            print()
        try:
            print("CASE 3b")
            ratio = (B**1.5*DX0**2*(X0 - X1)**2*((A**1.5*(U0 - U1)*(7*DX1*(U0 - U1) + (3*sqrt(A)*sqrt(B) - 4*original_DDX1)*(X0 - X1))*(DX1*(7*U0**2 - 14*U0*U1 + 7*U1**2 + 16*X0 - 16*X1) + 3*sqrt(A)*sqrt(B)*(U0 - U1)*(X0 - X1)))/(B**1.5*DX0**2*(X0 - X1)**2) - (12*(DX0*(U0 - U1)*(7*U0**2 - 14*U0*U1 + 7*U1**2 + 32*(X0 - X1)) + DX1*(U0 - U1)*(7*U0**2 - 14*U0*U1 + 7*U1**2 + 32*(X0 - X1)) + 2*(X0 - X1)*(3*sqrt(A)*sqrt(B)*(U0 - U1)**2 + 40*(-X0 + X1))))/(sqrt(A)*sqrt(B)*(X0 - X1)**2) + 2*sqrt(2)*sqrt((3*(U0 - U1)**2*(7*DX0*(U0 - U1) + 7*DX1*(U0 - U1) + 2*(3*sqrt(A)*sqrt(B) + 2*original_DDX0 - 2*original_DDX1)*(X0 - X1))*(A*DX1*(7*U0**2 - 14*U0*U1 + 7*U1**2 + 16*X0 - 16*X1) + 3*A**1.5*sqrt(B)*(U0 - U1)*(X0 - X1))**2 - 16*A**2.5*sqrt(B)*(DX1*(7*U0**2 - 14*U0*U1 + 7*U1**2 + 16*X0 - 16*X1) + 3*sqrt(A)*sqrt(B)*(U0 - U1)*(X0 - X1))**2*(X0 - X1)**2 - 3*A**2*(U0 - U1)*(7*DX1*(U0 - U1) + (3*sqrt(A)*sqrt(B) - 4*original_DDX1)*(X0 - X1))*(DX1*(7*U0**2 - 14*U0*U1 + 7*U1**2 + 16*X0 - 16*X1) + 3*sqrt(A)*sqrt(B)*(U0 - U1)*(X0 - X1))*(DX0*(U0 - U1)*(7*U0**2 - 14*U0*U1 + 7*U1**2 + 32*(X0 - X1)) + DX1*(U0 - U1)*(7*U0**2 - 14*U0*U1 + 7*U1**2 + 32*(X0 - X1)) + 2*(X0 - X1)*(3*sqrt(A)*sqrt(B)*(U0 - U1)**2 + 40*(-X0 + X1))) + 18*B*DX0**2*(DX0*(U0 - U1)*(7*U0**2 - 14*U0*U1 + 7*U1**2 + 32*(X0 - X1)) + DX1*(U0 - U1)*(7*U0**2 - 14*U0*U1 + 7*U1**2 + 32*(X0 - X1)) + 2*(X0 - X1)*(3*sqrt(A)*sqrt(B)*(U0 - U1)**2 + 40*(-X0 + X1)))**2)/(A*B**2*DX0**2*(X0 - X1)**4))))/(A**1.5*(DX1*(7*U0**2 - 14*U0*U1 + 7*U1**2 + 16*X0 - 16*X1) + 3*sqrt(A)*sqrt(B)*(U0 - U1)*(X0 - X1))**2)
            print_abg()
        except:
            print()
            print("Failed..")
            print()
        # CASE 4:  (gamma == bound) and (beta > 6)
        try:
            print("CASE 4a  (gamma == bound) and (beta > 6)")
            ratio = (sqrt(B)*DX0**2*(X0 - X1)**2*((sqrt(A)*(U0 - U1)*(7*DX0*(U0 - U1) + (3*sqrt(A)*sqrt(B) + 4*original_DDX0)*(X0 - X1))*(DX0*(7*U0**2 - 14*U0*U1 + 7*U1**2 + 16*X0 - 16*X1) + 3*sqrt(A)*sqrt(B)*(U0 - U1)*(X0 - X1)))/(sqrt(B)*DX0**2*(X0 - X1)**2) - (12*(DX0*(U0 - U1)*(7*U0**2 - 14*U0*U1 + 7*U1**2 + 32*(X0 - X1)) + DX1*(U0 - U1)*(7*U0**2 - 14*U0*U1 + 7*U1**2 + 32*(X0 - X1)) + 2*(X0 - X1)*(3*sqrt(A)*sqrt(B)*(U0 - U1)**2 + 40*(-X0 + X1))))/(sqrt(A)*sqrt(B)*(X0 - X1)**2) - 2*sqrt(2)*sqrt((3*A*(U0 - U1)**2*(7*DX0*(U0 - U1) + 7*DX1*(U0 - U1) + 2*(3*sqrt(A)*sqrt(B) + 2*original_DDX0 - 2*original_DDX1)*(X0 - X1))*(DX0*(7*U0**2 - 14*U0*U1 + 7*U1**2 + 16*X0 - 16*X1) + 3*sqrt(A)*sqrt(B)*(U0 - U1)*(X0 - X1))**2 - 16*A**1.5*sqrt(B)*(DX0*(7*U0**2 - 14*U0*U1 + 7*U1**2 + 16*X0 - 16*X1) + 3*sqrt(A)*sqrt(B)*(U0 - U1)*(X0 - X1))**2*(X0 - X1)**2 - 3*A*(U0 - U1)*(7*DX0*(U0 - U1) + (3*sqrt(A)*sqrt(B) + 4*original_DDX0)*(X0 - X1))*(DX0*(7*U0**2 - 14*U0*U1 + 7*U1**2 + 16*X0 - 16*X1) + 3*sqrt(A)*sqrt(B)*(U0 - U1)*(X0 - X1))*(DX0*(U0 - U1)*(7*U0**2 - 14*U0*U1 + 7*U1**2 + 32*(X0 - X1)) + DX1*(U0 - U1)*(7*U0**2 - 14*U0*U1 + 7*U1**2 + 32*(X0 - X1)) + 2*(X0 - X1)*(3*sqrt(A)*sqrt(B)*(U0 - U1)**2 + 40*(-X0 + X1))) + 18*DX0**2*(DX0*(U0 - U1)*(7*U0**2 - 14*U0*U1 + 7*U1**2 + 32*(X0 - X1)) + DX1*(U0 - U1)*(7*U0**2 - 14*U0*U1 + 7*U1**2 + 32*(X0 - X1)) + 2*(X0 - X1)*(3*sqrt(A)*sqrt(B)*(U0 - U1)**2 + 40*(-X0 + X1)))**2)/(A*B*DX0**2*(X0 - X1)**4))))/(sqrt(A)*(DX0*(7*U0**2 - 14*U0*U1 + 7*U1**2 + 16*X0 - 16*X1) + 3*sqrt(A)*sqrt(B)*(U0 - U1)*(X0 - X1))**2)
            print_abg()
        except:
            print()
            print("Failed..")
            print()
        try:
            print("CASE 4b  (gamma == bound) and (beta > 6)")
            ratio = (sqrt(B)*DX0**2*(X0 - X1)**2*((sqrt(A)*(U0 - U1)*(7*DX0*(U0 - U1) + (3*sqrt(A)*sqrt(B) + 4*original_DDX0)*(X0 - X1))*(DX0*(7*U0**2 - 14*U0*U1 + 7*U1**2 + 16*X0 - 16*X1) + 3*sqrt(A)*sqrt(B)*(U0 - U1)*(X0 - X1)))/(sqrt(B)*DX0**2*(X0 - X1)**2) - (12*(DX0*(U0 - U1)*(7*U0**2 - 14*U0*U1 + 7*U1**2 + 32*(X0 - X1)) + DX1*(U0 - U1)*(7*U0**2 - 14*U0*U1 + 7*U1**2 + 32*(X0 - X1)) + 2*(X0 - X1)*(3*sqrt(A)*sqrt(B)*(U0 - U1)**2 + 40*(-X0 + X1))))/(sqrt(A)*sqrt(B)*(X0 - X1)**2) + 2*sqrt(2)*sqrt((3*A*(U0 - U1)**2*(7*DX0*(U0 - U1) + 7*DX1*(U0 - U1) + 2*(3*sqrt(A)*sqrt(B) + 2*original_DDX0 - 2*original_DDX1)*(X0 - X1))*(DX0*(7*U0**2 - 14*U0*U1 + 7*U1**2 + 16*X0 - 16*X1) + 3*sqrt(A)*sqrt(B)*(U0 - U1)*(X0 - X1))**2 - 16*A**1.5*sqrt(B)*(DX0*(7*U0**2 - 14*U0*U1 + 7*U1**2 + 16*X0 - 16*X1) + 3*sqrt(A)*sqrt(B)*(U0 - U1)*(X0 - X1))**2*(X0 - X1)**2 - 3*A*(U0 - U1)*(7*DX0*(U0 - U1) + (3*sqrt(A)*sqrt(B) + 4*original_DDX0)*(X0 - X1))*(DX0*(7*U0**2 - 14*U0*U1 + 7*U1**2 + 16*X0 - 16*X1) + 3*sqrt(A)*sqrt(B)*(U0 - U1)*(X0 - X1))*(DX0*(U0 - U1)*(7*U0**2 - 14*U0*U1 + 7*U1**2 + 32*(X0 - X1)) + DX1*(U0 - U1)*(7*U0**2 - 14*U0*U1 + 7*U1**2 + 32*(X0 - X1)) + 2*(X0 - X1)*(3*sqrt(A)*sqrt(B)*(U0 - U1)**2 + 40*(-X0 + X1))) + 18*DX0**2*(DX0*(U0 - U1)*(7*U0**2 - 14*U0*U1 + 7*U1**2 + 32*(X0 - X1)) + DX1*(U0 - U1)*(7*U0**2 - 14*U0*U1 + 7*U1**2 + 32*(X0 - X1)) + 2*(X0 - X1)*(3*sqrt(A)*sqrt(B)*(U0 - U1)**2 + 40*(-X0 + X1)))**2)/(A*B*DX0**2*(X0 - X1)**4))))/(sqrt(A)*(DX0*(7*U0**2 - 14*U0*U1 + 7*U1**2 + 16*X0 - 16*X1) + 3*sqrt(A)*sqrt(B)*(U0 - U1)*(X0 - X1))**2)
            print_abg()
        except:
            print()
            print("Failed..")
            print()

        print()

        # low_bound,     upp_bound     = 0,    1
        # # Continue dividing the interval in 2 to find the smallest
        # # "upper bound" (amount target) that satisfies monotonicity.
        # while ((upp_bound - low_bound) > accuracy):
        #     to_target = upp_bound / 2  +  low_bound / 2
        #     # If we found the limit of floating point numbers, break.
        #     if ((to_target == upp_bound) or (to_target == low_bound)): break
        #     # Recompute DDX0 and DDX1 based on the to_target.
        #     DDX0 = (1-to_target) * original_DDX0 + to_target * target_DDX0
        #     DDX1 = (1-to_target) * original_DDX1 + to_target * target_DDX1
        #     # Otherwise, proceed with a binary seaarch.
        #     if is_monotone(): upp_bound = to_target
        #     else:             low_bound = to_target
        # # Store the smallest amount "target" for DDX0 and DDX1 that is monotone.
        # DDX0 = (1-upp_bound) * original_DDX0 + upp_bound * target_DDX0
        # DDX1 = (1-upp_bound) * original_DDX1 + upp_bound * target_DDX1

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


if __name__ == "__main__":
    # --------------------------------------------------------------------
    #               TEST CASES
    # 
    # 0, 4, (None)     -> Good
    # 1, 4, (None)     -> Bad (not monotone) (now good)
    # 0, 13, (None)    -> Bad (not monotone after first pass) (now good)
    # 0, 30, (-3,None) -> Bad (far right still not monotone after passes) (now good)
    # 0, 100, (-12,-7) -> Bad (fixing previous breaks second previous)
    # 
    # --------------------------------------------------------------------

    SEED = 0
    NODES = 100
    SUBSET = slice(-12,-7) # slice(None) 

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

    # Generate a plot to see what it all looks like.
    from util.plot import Plot
    p = Plot()
    p.add("Points", x, y)

    kwargs = dict(plot_points=1000)

    p.add_func("continuity 0", fit(x,y,continuity=0, non_decreasing=True),
               [min(x), max(x)], group='0', **kwargs)
    p.add_func("c0 d1", fit(x,y,continuity=0, non_decreasing=True).derivative(),
               [min(x), max(x)], dash="dash", color=p.color(p.color_num,alpha=.5), 
               group='0', **kwargs)

    p.add_func("continuity 1", fit(x,y,continuity=1,
                                   non_decreasing=False, ends=1, mids=0), 
               [min(x), max(x)], group='1', **kwargs)
    p.add_func("c1 d1", fit(x,y,continuity=1, non_decreasing=True).derivative(), 
               [min(x), max(x)], dash="dash", color=p.color(p.color_num,alpha=.5), 
               group='1', **kwargs)

    f = monotone_cubic_spline(x,y)
    p.add_func("monotone c1", f, [min(x), max(x)], group='1m', **kwargs)
    p.add_func("monotone c1 d1", f.derivative(), [min(x), max(x)], 
               dash="dash", color=p.color(p.color_num,alpha=.5), group='1m', **kwargs)

    p.add_func("continuity 2", fit(x,y,continuity=2, non_decreasing=True), 
               [min(x), max(x)], group='2', **kwargs)
    p.add_func("c2 d1", fit(x,y,continuity=2, non_decreasing=True).derivative(), 
               [min(x), max(x)], dash="dash", color=p.color(p.color_num,alpha=.5), group='2', **kwargs)

    f = monotone_quintic_spline(x,y, fix_previous=False)
    p.add_func("monotone c2 (no fix)", f, [min(x), max(x)], group='2m', color=p.color(6), **kwargs)
    p.add_func("monotone c2d1 (no fix)", f.derivative(), [min(x), max(x)], 
               dash="dash", color=p.color(6,alpha=.5), group='2m', **kwargs)

    f = monotone_quintic_spline(x,y)
    p.add_func("monotone c2", f, [min(x), max(x)], group='2mf', color=p.color(7), **kwargs)
    p.add_func("monotone c2d1", f.derivative(), [min(x), max(x)], 
               dash="dash", color=p.color(7,alpha=.5), group='2mf', **kwargs)

    p.show(file_name="monotone_quintic_interpolating_spline.html")
