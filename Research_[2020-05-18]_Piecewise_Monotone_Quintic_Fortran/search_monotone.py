from polynomial import Spline, polynomial

# Convert a value or list of values into a list of Fraction objects.
def make_exact(value):
    from fraction import Fraction
    try:    return list(map(make_exact, value))
    except: return Fraction(value)

# Compute a monotone quintic spline that interpolates given data.
def monotone_quintic_spline(x, y, exact=True, accuracy=2**(-26),
                            verbose=False, local_fits=None,
                            free=False, monotone=True):
    # Initialize dicationary of local fits.
    if local_fits is None: local_fits = {}
    # Make the values exact.
    if exact:
        x = make_exact(x)
        y = make_exact(y)
    # Identify the local extreme points.
    extremes = {i for i in range(1, len(x)-1)
                      if ((y[i]-y[i-1]) * (y[i+1]-y[i]) < 0)}
    flats = {i for i in range(1, len(x)-1)
                   if ((y[i]-y[i-1]) * (y[i+1]-y[i]) == 0)}
    # Compute "ideal" values everywhere.
    values = [[yi] for yi in y]
    # Compute first and second derivatives.
    for i in range(len(x)): values[i] += estimate_derivative(
            x, y, i, extremes, flats, free=free)
    # Store local fits.
    for i in range(len(x)):
        local_fits[i] = (
            polynomial([x[i]], [y[i]], dx0=values[i][1], ddx0=values[i][2]),
            [x[max(0,i-1)], x[min(i+1,len(x)-1)]] )
    # Find values nearest to the ideal values that are monotone.
    if monotone:
        # Store the current "step size" that was used for all values.
        best_values = [vi.copy() for vi in values]
        # Make step size and values exact.
        if exact: best_values, values = make_exact(best_values), make_exact(values)
        # Track which ones have been shrank.
        did_shrink = {}
        # Walk as close to the ideal values as possible.
        step = 0
        step_size = 1
        # Identify those intervals that need to be made monotone.
        if (not free): monotone_intervals = list(range(len(x)-1))
        else: monotone_intervals = [i for i in range(len(x)-1) if (
                (values[i][1] * (values[i+1][0] - values[i][0]) >= 0) and
                (values[i+1][1] * (values[i+1][0] - values[i][0]) >= 0))]
        while (step_size > accuracy) or any(
                not is_monotone(x[i], x[i+1], *values[i], *values[i+1])
                for i in monotone_intervals):
            # Set the step size for this iteration.
            step_size = max(step_size / 2, accuracy)
            step += 1
            # Find the values that need to be shrunk.
            to_shrink = {}
            # Step any intervals that are not monotone backwards by 'step size'.
            for i in monotone_intervals:
                if not is_monotone(x[i], x[i+1], *values[i], *values[i+1]):
                    to_shrink[i] = True
                    to_shrink[i+1] = True
            # Shrink those values that need to be shrunk.
            for i in sorted(to_shrink):
                # Shrink the first derivative (bounded by 0).
                values[i][1] = values[i][1] - step_size * best_values[i][1]
                if (best_values[i][1] * values[i][1] < 0): values[i][1] = 0
                # Shrink the second derivative (bounded by 0).
                values[i][2] = values[i][2] - step_size * best_values[i][2]
                if (best_values[i][2] * values[i][2] < 0): values[i][2] = 0
                # Record that this value has been shrunk.
                did_shrink[i] = True
            # Grow any values that were shrunk, but not in this iteration.
            for i in did_shrink:
                if i in to_shrink: continue
                if (step_size <= accuracy): continue
                # Grow the first derivative (bounded by ideal value).
                values[i][1] = values[i][1] + step_size * best_values[i][1]
                sign = (-1) ** (best_values[i][1] < 0)
                if (sign * (best_values[i][1] - values[i][1]) < 0):
                    values[i][1] = best_values[i][1]
                # Grow the second derivative (bounded by ideal value).
                values[i][2] = values[i][2] + step_size * best_values[i][2]
                sign = (-1) ** (best_values[i][2] < 0)
                if (sign * (best_values[i][2] - values[i][2]) < 0):
                    values[i][2] = best_values[i][2]
    # Return the monotone quintic spline.
    return Spline(x, values)
    

# Given "i" the index at which the first derivative should be
# estimated, construct local quadratic fits and pick the slope of the
# one with the lowest curvature.
def estimate_derivative(x, y, i, extremes, flats, free=False):
    # If this is a local flat, estimate 0.
    if (i in flats): return [0, 0]
    # If this is a local maximum, force first derivative to zero.
    elif (i in extremes):
        direction = 0.0
        functions = []
        # Compute the left function (interpolates zero slope).
        if (i > 0): functions.append( polynomial(
                [x[i-1], x[i]], [y[i-1],y[i]], dx1=0) )
        # Compute the right function (interpolates zero slope).
        if (i < len(x)-1): functions.append( polynomial(
                [x[i], x[i+1]], [y[i],y[i+1]], dx0=0) )
    else:
        # Construct quadratic interpolants over a sliding window.
        functions = []
        # Add the appropriate left interpolant (i-2, i-1, i).
        if (i > 0):
            if (i-1 in extremes) or (i-1 in flats):
                f = polynomial([x[i-1],x[i]], [y[i-1],y[i]], dx0=0)
                functions.append(f)
            elif (i-2 >= 0):
                if (y[i-2] == y[i-1]):
                    f = polynomial([x[i-1],x[i]], [y[i-1],y[i]], dx0=0)
                    functions.append(f)
                else:
                    f = polynomial(x[i-2:i+1], y[i-2:i+1])
                    functions.append(f)
        # Add the appropriate center interpolant (i-1, i, i+1).
        if ((i-1 >= 0) and (i+1 < len(x)) and not
            (((i-1 in flats) or (i-1 in extremes)) and
             ((i+1 in flats) or (i+1 in extremes)))):
            f = polynomial(x[i-1:i+2], y[i-1:i+2])
            functions.append(f)
        # Add the appropriate right interpolant (i, i+1, i+2).
        if (i + 1 < len(x)):
            if (i+1 in extremes) or (i+1 in flats):
                f = polynomial([x[i],x[i+1]], [y[i],y[i+1]], dx1=0)
                functions.append(f)
            elif (i+2 < len(x)):
                f = polynomial(x[i:i+3], y[i:i+3])
                functions.append(f)
        # Set the direction.
        if (i > 0):
            if (y[i-1] < y[i]): direction = 1.0
            elif (y[i-1] > y[i]): direction = -1.0
        elif (i+1 < len(x)):
            if (y[i] < y[i+1]): direction = 1.0
            elif (y[i] > y[i+1]): direction = -1.0
    # Sort the functions by their curvature.
    derivatives = []
    for f in functions:
        df = f.derivative()
        dxi = df(x[i])
        # Skip estimates with derivatives that are pointing the wrong way.
        if (dxi * direction < 0): continue
        ddf = df.derivative()
        ddxi = ddf(x[i])
        derivatives += [(dxi, ddxi)]
    # If there were no viable options, return 0's.
    if (len(derivatives) == 0): return [0, 0]
    # Sort the derivatives by magnitude of curvature, then by
    # magnitude of first derivative when curvatures are equal.
    derivatives.sort(key=lambda d: abs(d[1]))
    return derivatives[0]


# Function for computing tight condition on monotonicity. Returns True
# if an interval is monotone, False otherwise.
def is_monotone(U0, U1, X0, DX0, DDX0, X1, DX1, DDX1):
    # Flip the sign to only consider the monotone increasing case.
    sign = (-1) ** int(X1 < X0)
    X0 *= sign
    X1 *= sign
    DX0 *= sign
    DX1 *= sign
    DDX0 *= sign
    DDX1 *= sign
    # Make sure the slopes point in the right direction.
    if ((X1 - X0) * DX0) < 0: return False
    if ((X1 - X0) * DX1) < 0: return False
    # Compute A and B.
    A = (U1 - U0) * DX0 / (X1 - X0)
    B = (U1 - U0) * DX1 / (X1 - X0)
    # Simplified cubic monotone case.
    if (A*B <= 0):
        alpha = (4*DX1 + DDX1*(U0-U1)) * (U0-U1) / (X0-X1)
        beta = 30 + ((-24*(DX0+DX1) + 3*(DDX0-DDX1)*(U0-U1))*(U0-U1)) / (2 * (X0-X1))
        gamma = ((U0-U1) * (4*DX0 + DDX0*(U1-U0))) / (X0-X1)
        delta = DX0 * (U0-U1) / (X0-X1)
        return (alpha >= 0) and (delta >= 0) and \
            (beta >= alpha - (4*alpha*delta)**(1/2)) and \
            (gamma >= delta - (4*alpha*delta)**(1/2))
    # Full quintic monotone case.
    tau_1 = 24 + 2*(A*B)**(1/2) - 3*(A+B)
    if (tau_1 < 0): return False
    # Compute DDX0 and DDX1 that satisfy monotonicity by scaling (C,D)
    # down until montonicity is achieved (using binary search).
    alpha_constant   = 4 * (B**(1/4) / A**(1/4))
    alpha_multiplier = ((U0-U1) / DX1) * B**(1/4) / A**(1/4)
    gamma_constant   = 4 * (DX0 / DX1) * (B**(3/4) / A**(3/4))
    gamma_multiplier = ((U1-U0) / DX1) * (B**(3/4) / A**(3/4))
    beta_constant    = (12 * (DX0+DX1) * (U1-U0) + 30 * (X0-X1)) / ((X0-X1) * A**(1/2) * B**(1/2))
    beta_multiplier  = (3 * (U0-U1)**2) / (2 * (X0-X1) * A**(1/2) * B**(1/2))
    # Compute the monotonicity condition.
    a = alpha_constant + alpha_multiplier * DDX1
    g = gamma_constant + gamma_multiplier * DDX0
    b = beta_constant  + beta_multiplier  * (DDX0 - DDX1)
    if b <= 6: bound = - (b + 2) / 2
    else:      bound = -2 * (b - 2)**(1/2)
    return (a > bound) and (g > bound)
