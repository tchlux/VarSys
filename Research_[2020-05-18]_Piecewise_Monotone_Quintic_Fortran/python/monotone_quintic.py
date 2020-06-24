from fraction import Fraction
from polynomial import Spline, polynomial

# Convert a value or list of values into a list of Fraction objects.
def make_exact(value):
    try:    return list(map(make_exact, value))
    except: return Fraction(value)

# Compute a monotone quintic spline that interpolates given data.
def monotone_quintic_spline(x, y, accuracy=2**(-26), verbose=False, monotone=True):
    # Make the values exact.
    x, y = make_exact(x), make_exact(y)
    # Identify the local extreme points.
    extremes = {i for i in range(1, len(x)-1)
                      if ((y[i]-y[i-1]) * (y[i+1]-y[i]) < 0)}
    flats = {i for i in range(1, len(x)-1)
                   if ((y[i]-y[i-1]) * (y[i+1]-y[i]) == 0)}
    # Compute "ideal" values everywhere.
    values = [[yi] for yi in y]
    # Compute first and second derivatives.
    for i in range(len(x)): values[i] += quadratic_facet(
            x, y, i, extremes, flats)
    # Find values nearest to the ideal values that are monotone.
    if monotone:
        # Store the current "step size" that was used for all values.
        best_values = [vi.copy() for vi in values]
        # Make step size and values exact.
        best_values, values = make_exact(best_values), make_exact(values)
        # Walk as close to the ideal values as possible.
        step_size = 1
        # Track three queues:
        to_check = {}   # Intervals that need to be checked.
        to_shrink = {}  # Intervals that need to be corrected.
        to_grow = {}    # Intervals that were previously corrected.
        # Find the values that need to be shrunk.
        for i in range(len(x)-1):
            if not is_monotone(x[i], x[i+1], *values[i], *values[i+1]):
                to_shrink[i] = True
                to_shrink[i+1] = True
        # Record whether or not the minimum step has been reached.
        reached_min_step = False
        while (not reached_min_step) or (len(to_shrink) > 0):
            # Set the step size for this iteration.
            if (not reached_min_step):
                step_size /= 2
                if (step_size < accuracy):
                    step_size = accuracy
                    reached_min_step = True
                    to_grow = {}
            # If the minimum step size was reached, start increasing again.
            else: step_size = step_size + step_size / 2
            # Grow any values that were shrunk, but not in this iteration.
            for i in to_grow:
                if i in to_shrink: continue
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
                # Record that this interval needs to be checked.
                if (i > 0):          to_check[i-1] = True
                if (i < (len(x)-1)): to_check[i] = True
            # Shrink those values that need to be shrunk.
            for i in to_shrink:
                # Shrink the first derivative (bounded by 0).
                values[i][1] = values[i][1] - step_size * best_values[i][1]
                if (best_values[i][1] * values[i][1] < 0): values[i][1] = 0
                # Shrink the second derivative (bounded by 0).
                values[i][2] = values[i][2] - step_size * best_values[i][2]
                if (best_values[i][2] * values[i][2] < 0): values[i][2] = 0
                # Record that this value was shrunk and should be checked.
                if (not reached_min_step): to_grow[i] = True
                if (i > 0):          to_check[i-1] = True
                if (i < (len(x)-1)): to_check[i] = True
            # Find the values that need to be shrunk.
            to_shrink = {}
            # Step any intervals that are not monotone backwards by 'step size'.
            for i in to_check:
                if not is_monotone(x[i], x[i+1], *values[i], *values[i+1]):
                    to_shrink[i] = True
                    to_shrink[i+1] = True
            to_check = {}
    # Return the monotone quintic spline.
    return Spline(x, values)
    

# Given "i" the index at which the first derivative should be
# estimated, construct local quadratic fits and pick the slope of the
# one with the lowest curvature.
def quadratic_facet(x, y, i, extremes, flats):
    # If this is a local flat, estimate 0.
    if (i in flats): return [0, 0]
    # If this is a local maximum, force first derivative to zero, 
    # and assume that it is not an endpoint.
    elif (i in extremes):
        direction = 0
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
        direction = 0
        if (i > 0):
            if (y[i-1] < y[i]):   direction =  1
            elif (y[i-1] > y[i]): direction = -1
        elif (i+1 < len(x)):
            if (y[i] < y[i+1]):   direction =  1
            elif (y[i] > y[i+1]): direction = -1
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
    # Sort the derivatives by magnitude of curvature, then return the 
    # values from the function with the least curvature.
    derivatives.sort(key=lambda d: abs(d[1]))
    return derivatives[0]


# Function for computing tight condition on monotonicity. Returns True
# if an interval is monotone, False otherwise.
def is_monotone(U0, U1, F0, DF0, DDF0, F1, DF1, DDF1):
    # Flip the sign to only consider the monotone increasing case.
    if (F1 < F0):
        F0,DF0,DDF0 = -F0, -DF0, -DDF0
        F1,DF1,DDF1 = -F1, -DF1, -DDF1
    # Make sure the slopes point in the right direction.
    if ((F1 - F0) * DF0) < 0: return False
    if ((F1 - F0) * DF1) < 0: return False
    # Compute A and B.
    A = (U1 - U0) * DF0 / (F1 - F0)
    B = (U1 - U0) * DF1 / (F1 - F0)
    # Simplified cubic monotone case.
    if (A*B <= 0):
        alpha = (4*DF1 + DDF1*(U0-U1)) * (U0-U1) / (F0-F1)
        beta = 30 + ((-24*(DF0+DF1) + 3*(DDF0-DDF1)*(U0-U1))*(U0-U1)) / (2 * (F0-F1))
        gamma = ((U0-U1) * (4*DF0 + DDF0*(U1-U0))) / (F0-F1)
        delta = DF0 * (U0-U1) / (F0-F1)
        return (alpha >= 0) and (delta >= 0) and \
            (beta >= alpha - (4*alpha*delta)**(1/2)) and \
            (gamma >= delta - (4*alpha*delta)**(1/2))
    # Full quintic monotone case.
    tau_1 = 24 + 2*(A*B)**(1/2) - 3*(A+B)
    if (tau_1 < 0): return False
    alpha = (4*DF1 + DDF1*(U0-U1)) / (DF0 * DF1**3)**(1/4)
    gamma = (4*DF0 + DDF0*(U1-U0)) / (DF0**3 * DF1)**(1/4)
    beta = (60*(F1-F0)/(U1-U0) + 3*((DDF1-DDF0)*(U1-U0) - 8*(DF0+DF1))) / (2*(DF0*DF1)**(1/2))
    if beta <= 6: bound = -(beta + 2) / 2
    else:         bound = -2 * (beta - 2)**(1/2)
    return (alpha > bound) and (gamma > bound)
