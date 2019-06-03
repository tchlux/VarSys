
# Given unique "x" values and associated "y" values (of same length),
# construct an interpolating polynomial with the Newton divided
# difference method. Return that polynomial as a function.
def polynomial(x, y):
    # Sort the data by "x" value.
    indices = sorted(range(len(x)), key=lambda i: x[i])
    x = [x[i] for i in indices]
    y = [y[i] for i in indices]
    # Compute the divided difference table.
    dd_values = [y]
    for d in range(1, len(x)):
        slopes = []
        for i in range(len(dd_values[-1])-1):
            try:    dd = (dd_values[-1][i+1] - dd_values[-1][i]) / (x[i+d] - x[i])
            except: dd = 0
            slopes.append( dd )
        dd_values.append( slopes )
    # Get the divided difference (polynomial coefficients) in reverse
    # order so that the most nested value is first.
    coefs = [row[0] for row in reversed(dd_values)]
    points = list(reversed(x))
    # Generate a polynomial function (with numerically stable
    #  evaluation) from the set of zeros and coefficients.
    def f(z, order=len(coefs)):
        total = coefs[0]
        for d in range(1,order):
            total = coefs[d] + (z - points[d]) * total
        return total
    f.points = points
    f.coefs = coefs
    # Return the interpolating polynomial.
    return f

# Given a (left value, left d1, ...), (right value, right d1, ...)
# pair of tuples, return the lowest order polynomial necessary to
# exactly match those values and derivatives at 0 on the left and 1
# on the right (any range can be mapped into this).
def polynomial_piece(left, right, interval=(0,1)):
    # Make sure both are lists.
    left  = list(left)
    right = list(right)
    # Fill values by matching them on both sides of interval (reducing order).
    for i in range(len(left) - len(right)):
        right.append( left[len(right)] )
    for i in range(len(right) - len(left)):
        left.append( right[len(left)] )
    # First match the function value, then compute the coefficients
    # for all of the higher order terms in the polynomial.
    coefs = list(left)
    # Compute the divided difference up until we flatten out the
    # highest provided derivative between the two sides.
    dds = [ [left[0], right[0]], [right[0] - left[0]] ]
    for i in range(len(left)-2):
        new_vals = [ dds[-1][0] - left[i+1] ]
        for (d_left, d_right) in zip(dds[-1][:-1], dds[-1][1:]):
            new_vals.append( d_right - d_left )
        new_vals.append( right[i+1] - dds[-1][-1] )
        dds.append( new_vals )
    # Now the last row of the dd table should be level with "left" and
    # "right", we can complete the rest of the table in the normal way.
    row = [left[-1]] + dds[-1] + [right[-1]]
    while (len(row) > 1):
        row = [(r-l) for (l,r) in zip(row[:-1], row[1:])]
        coefs.append(row[0])
    # Remove the last additional element for linear functions.
    if (len(coefs) == 3): coefs.pop(-1)
    # Reverse the coefficients to go from highest order (most nested)
    # to lowest order, making evaluation slightly more clean.
    points = [1]*len(left) + [0]*len(left)
    coefs = list(reversed(coefs))
    shift = interval[0]
    scale = interval[1] - interval[0]
    # Generate a polynomial function (with numerically stable
    #  evaluation) from the set of zeros and coefficients.
    def f(z):
        # Rescale "z" to the interpolation interval for this polynomial.
        z = (z - shift) / scale
        # Compute the value (in a numerically stable way).
        total = coefs[0]
        for d in range(1,len(coefs)):
            total = coefs[d] + (z-points[d]) * total
        return total
    f.points = points
    f.coefs = coefs
    # Return the polynomial function.
    return f

# Approximate a polynomial piece using a very small number offset and
# a full Newton interpolating polynomial.
def approx_piece(left, right, step=.001):
    # Make sure both are lists.
    left  = list(left)
    right = list(right)
    # Fill values by matching them on both sides of interval (reducing order).
    for i in range(len(left) - len(right)):
        right.append( left[len(right)] )
    for i in range(len(right) - len(left)):
        left.append( right[len(left)] )
    # Generate a set of points to approximate the Newton polynomial source.
    x = []
    y = []
    left_func = lambda x: sum(x**i * v for (i,v) in enumerate(left))
    right_func = lambda x: sum((x-1)**i * v for (i,v) in enumerate(right))
    for i in range(len(left)):
        x.append( step*i)
        y.append( left_func(x[-1]) )
    for i in range(len(right)-1,-1,-1):
        x.append( 1 - step*i)
        y.append( right_func(x[-1]) )
    return polynomial(x,y)


if __name__ == "__main__":
    TEST_NEWTON_INTERPOLANT = True
    if TEST_NEWTON_INTERPOLANT:
        SMALL = 1.4901161193847656*10**(-8) 
        # ^^ SQRT(EPSILON(REAL(1.0)))
        x_vals = [0,1,2,3,4,5]
        y_vals = [1,2,1,2,1,10]
        f = polynomial(x_vals, y_vals)
        for (x,y) in zip(x_vals,y_vals):
            try:    assert( abs(y - f(x)) < SMALL )
            except:
                string =  "\n\nFailed test.\n"
                string += f" x:    {x}\n"
                string += f" y:    {y}\n"
                string += f" f({x}): {f(x)}"
                class FailedTest(Exception): pass
                raise(FailedTest(string))

    PLOT = False
    TEST_POLY_APPROX = True
    if TEST_POLY_APPROX:
        if PLOT:
            from util.plot import Plot
            p = Plot("Polynomial Pieces")
        # Pick the (value, d1, ...) pairs for tests.
        left_rights = [
            ([0], [0]),
            ([0], [1]),
            ([0,1], [0,-1]),
            ([0,2], [0,-2]),
            ([0,1], [1,0]),
            ([1,0], [0,-1]),
            ([0,1], [0,1]),
            ([0,1,0], [0,-1,1]),
            ([10,7,3], [19,3]),
        ]
        # Plot a bunch of sample functions.
        for (left, right) in left_rights:
            f = approx_piece( left, right )
            approx_coefs = list(map(round,f.coefs))[::-1]
            approx_points = list(map(round,f.points))[::-1]
            f = polynomial_piece( left, right )
            exact_coefs = list(map(round,f.coefs))[::-1]
            exact_points = list(map(round,f.points))[::-1]
            # Verify that the approximation points are correct.
            try: assert(approx_coefs == exact_coefs)
            except:
                string =  "\n\nFailed test.\n"
                string += f" Approximate coefficients:   {approx_coefs}\n"
                string += f" Exact coefficients (wrong): {exact_coefs}\n"
                class FailedTest(Exception): pass
                raise(FailedTest(string))
            # Verify that the interpolation points are correct.
            try: assert(approx_points == exact_points)
            except:
                string =  "\n\nFailed test.\n"
                string += f" Approximate points:   {approx_points}\n"
                string += f" Exact points (wrong): {exact_points}\n"
                class FailedTest(Exception): pass
                raise(FailedTest(string))
            # Create a plot of the functions if a demo is desired.
            if PLOT: p.add_func(f"{left}  {right}", f, [-.1, 1.1])
        if PLOT: p.show(file_name="piecewise_polynomial.html")


    TEST_POLY_PIECE = False
    if TEST_POLY_PIECE:
        from util.plot import Plot
        p = Plot("Newton Interpolation")
        p.add("Points", x_vals, y_vals)
        p.add_func("Interpolant", f, [min(x_vals), max(x_vals)])
        p.show(file_name="divided_diff.html", show=False)

        p = Plot("Polynomial Piece")

        left, right = [0], [0]
        f = polynomial_piece( left, right )
        p.add_func("Line from (0) to (0)", f, [-.1, 1.1])

        left, right = [0], [1]
        f = polynomial_piece( left, right )
        p.add_func("Line from (0) to (1)", f, [-.1, 1.1])

        left, right = [0,1], [1,0]
        f = polynomial_piece( left, right )
        p.add_func("Quad from (0,1) to (1,0)", f, [-.1, 1.1])

        left, right = [1,0], [0,-1]
        f = polynomial_piece( left, right )
        p.add_func("Quad from (1,0) to (0,-1)", f, [-.1, 1.1])

        left, right = [0,1,0], [0,-1,0]
        f = polynomial_piece( left, right )
        p.add_func("Quad from (0,1) to (0,-1)", f, [-.1, 1.1])

        p.show(file_name="divided_diff.html", append=True)
