# Fortran code to ...
#   create a Hermite interpolating polynomial out of B-splines.
#   create a Hermite interpolating spline out of B-splines.
#   compute the L2 difference between two splines.
#   test for monotonnicity of a spline.
#   identify the furthest monotone spline from 0 in direction.



# Given a list of knots and a list of values, find a good monotone
# interpolant for the data.
def monotone_fit(knots, values, continuity=2, ends=2, mids=2):
    print("Finding monotone fit..")
    from polynomial import Spline, fit
    # Construct an initial fit that is twice continuous.
    f = fit(x, y, continuity=continuity, non_decreasing=True, ends=ends, mids=mids)
    values = f.values
    values = np.array([list(map(float,vals)) for vals in values])
    values = nearest_monotone(knots, values)
    print(" done.")
    return Spline(knots, [list(map(Fraction,vals)) for vals in values])


# Given a list of knots and a list of values for the function and
# derivatives, walk along the line between all zero values and the
# provided values until the point where monotonicity is broken.
def nearest_monotone(knots, values):
    # Convert all inputs into numpy arrays.
    import numpy as np
    knots = np.array(knots)
    values = np.array(values)
    # Construct an array where all derivatives are assigned to zero.
    zero = np.zeros(values.shape)
    zero[:,0] = values[:,0]
    # Construct a test function that searches the line between the
    # zero values and the set of values provided.
    def test_on_line_function(ratio):
        return is_monotone_spline(knots, (1-ratio)*zero + ratio*values)
    # Find the nearest point on that line.
    ratio = binary_search(test_on_line_function, 0, 1)
    # Return the nearest monotone set of values.
    return (1-ratio)*zero + ratio*values


# Perform a binary search on a line over a boolean objective function.
# Find the point at which the the function switches outputs.
# 
# INPUTS:
#   f -- A boolean function of the form f(x) = (x <= z) for (l <= z < u).
#   l -- The lower bound for the range to search, f(l) = True.
#   u -- The upper bound for the range to search, f(u) = False.
# 
# OPTIONAL INPUTS:
#   accuracy -- The maximum allowable distance from the correct solution.
#   round    -- A rounding function to set the values of l and u.
# 
# RETURNS:
#   x -- A point within `accuracy` of the position where f transitions
#        from True to False.
#    
def binary_search(f, l, u, accuracy=2**(-26), round=lambda x: x):
    # If the top is valid, then return it.
    if f(u): return u
    # Verify the inputs provided are valid.
    assert(f(l) == True)
    assert(l <= u)
    assert(accuracy > 0)
    # Perform the binary search until we reach desired accuracy.
    while (abs(u - l) > accuracy):
        mid = round(l/2 + u/2)
        # If we've reached the limit of achievable accuracy, return.
        if (mid == l) or (mid == u): return mid
        # Transition the search window.
        if f(mid): l = mid
        else:      u = mid
    # Return the largest True locationn found for `f` on the intervale [l,u].
    return l

# Given a set of knots and values that define a hermite
# interpolating spline, return True if the spline is monotone.
def is_monotone_spline(knots, values, eval_points=1000):
    from polynomial import Spline
    # Construct the exact Hermite interpolating polynomial for this interval.
    f = Spline(knots, [list(map(Fraction,vals)) for vals in values])
    # Evaluate the interpolating spline at many points, see if it is monotone.
    evaluation_points = np.linspace(knots[0], knots[-1], 1000)
    evaluations = np.array(list(map(f, evaluation_points)))
    difference = evaluations[1:] - evaluations[:-1]
    signs = set(np.sign(difference)).intersection({1,-1})
    return len(signs) == 1




# 
if __name__ == "__main__":
    from polynomial import fit
    SEED = 0
    NODES = 13
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
    p.color_num += 1
    kwargs = {}
    # Cycle through all the different levels of continuity.
    for continuity in range(0,3+1):
        # Add the three functions to the plot.
        f = fit(x, y, continuity=continuity, non_decreasing=True)
        p.add_func(f"continuity {continuity}", f, interval,
                   group=continuity, dash="dot", **kwargs)
        mf = monotone_fit(x, y, continuity=continuity)
        p.add_func(f"monotone c{continuity}", mf, interval,
                   color=p.color(p.color_num), group=continuity, **kwargs)
        p.add_func(f"monotone c{continuity} derivative", mf.derivative(), interval, dash="dash",
                   color=p.color(p.color_num, alpha=.5), group=continuity, **kwargs)

    # Add the points last (so they are on top of everything).
    p.add("Points", list(map(float, x)), list(map(float, y)), color=p.color(0))
    # Generate the visual.
    p.show()
