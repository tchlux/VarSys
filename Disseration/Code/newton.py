
# Given unique "x" values and associated "y" values (of same length),
# construct an interpolating polynomial with the Newton divided
# difference method. Return that polynomial as a function.
def newton_interpolant(x, y):
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
    x = list(reversed(x))
    # Generate a polynomial function (with numerically stable
    #  evaluation) from the set of zeros and coefficients.
    def f(z, order=len(coefs)):
        total = coefs[0]
        for d in range(1,order):
            total = coefs[d] + (z-x[d]) * total
        return total
    f.points = x
    f.coefs = coefs
    # Return the interpolating polynomial.
    return f

# Given a (left value, left d1, ...), (right value, right d1, ...)
# pair of tuples, return the lowest order polynomial necessary to
# exactly match those values and derivatives at 0 on the left and 1
# on the right (any range can be mapped into this).
def polynomial_piece(left, right):
    # Make sure both are lists.
    left  = list(left)
    right = list(right)
    # Fill values by matching them on both sides of interval (reducing order).
    for i in range(len(left) - len(right)):
        right.append( left[len(right) + i] )
    for i in range(len(right) - len(left)):
        left.append( right[len(left) + i] )
    # First match the function value, then compute the coefficients
    # for all of the higher order terms in the polynomial.
    coefs = [left[0]]
    divisor = 1
    for i in range(len(left)):
        # Increment the divisor appropriately.
        if (i >= 2): divisor *= i
        # Get the left and right values.
        left_val = left[i]
        right_val = right[i]
        # Compute the coefficient.
        coefs.append( (right_val - left_val) / divisor )
    # Reverse the coefficients to go from highest order (most nested)
    # to lowest order, making evaluation slightly more clean.
    coefs = list(reversed(coefs))
    # Generate a polynomial function (with numerically stable
    #  evaluation) from the set of zeros and coefficients.
    def f(z):
        total = coefs[0]
        for d in range(1,len(coefs)):
            total = coefs[d] + (z-1) * total
        return total
    f.points = [1] * (len(coefs) - 1)
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
        right.append( left[len(right) + i] )
    for i in range(len(right) - len(left)):
        left.append( right[len(left) + i] )
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
    return newton_interpolant(x,y)


if __name__ == "__main__":
    TEST_NEWTON_INTERPOLANT = True
    if TEST_NEWTON_INTERPOLANT:
        SMALL = 1.4901161193847656*10**(-8) 
        # ^^ SQRT(EPSILON(REAL(1.0)))
        x_vals = [0,1,2,3,4,5]
        y_vals = [1,2,1,2,1,10]
        f = newton_interpolant(x_vals, y_vals)
        for (x,y) in zip(x_vals,y_vals):
            try:    assert( abs(y - f(x)) < SMALL )
            except:
                string =  "\n\nFailed test.\n"
                string += f" x:    {x}\n"
                string += f" y:    {y}\n"
                string += f" f({x}): {f(x)}"
                class FailedTest(Exception): pass
                raise(FailedTest(string))

    TEST_POLY_APPROX = True
    if TEST_POLY_APPROX:
        from util.plot import Plot
        p = Plot("Polynomial Pieces")

        print()
        print()
        print("------- Linear -------")
        left, right = [0], [0]
        print()
        print("left:    ",left)
        print("right:   ",right)
        f = approx_piece( left, right )
        print("f.coefs: ",list(map(round,f.coefs))[::-1])
        p.add_func(f"{left}  {right}", f, [-.1, 1.1])

        left, right = [0], [1]
        print()
        print("left:    ",left)
        print("right:   ",right)
        f = approx_piece( left, right )
        print("f.coefs: ",list(map(round,f.coefs))[::-1])
        p.add_func(f"{left}  {right}", f, [-.1, 1.1])

        print()
        print()
        print("------- Cubic -------")
        left, right = [0,1], [0,-1]
        print()
        print("left:    ",left)
        print("right:   ",right)
        f = approx_piece( left, right )
        print("f.coefs: ",list(map(round,f.coefs))[::-1])
        p.add_func(f"{left}  {right}", f, [-.1, 1.1])

        left, right = [0,2], [0,-2]
        print()
        print("left:    ",left)
        print("right:   ",right)
        f = approx_piece( left, right )
        print("f.coefs: ",list(map(round,f.coefs))[::-1])
        p.add_func(f"{left}  {right}", f, [-.1, 1.1])

        left, right = [0,1], [1,0]
        print()
        print("left:    ",left)
        print("right:   ",right)
        f = approx_piece( left, right )
        print("f.coefs: ",list(map(round,f.coefs))[::-1])
        p.add_func(f"{left}  {right}", f, [-.1, 1.1])

        left, right = [1,0], [0,-1]
        print()
        print("left:    ",left)
        print("right:   ",right)
        f = approx_piece( left, right )
        print("f.coefs: ",list(map(round,f.coefs))[::-1])
        p.add_func(f"{left}  {right}", f, [-.1, 1.1])

        left, right = [0,1], [0,1]
        print()
        print("left:    ",left)
        print("right:   ",right)
        f = approx_piece( left, right )
        print("f.coefs: ",list(map(round,f.coefs))[::-1])
        p.add_func(f"{left}  {right}", f, [-.1, 1.1])

        print()
        print()
        print("------- Quintic -------")
        left, right = [0,1,0], [0,-1,1]
        print()
        print("left:    ",left)
        print("right:   ",right)
        f = approx_piece( left, right )
        print("f.coefs: ",list(map(round,f.coefs))[::-1])
        p.add_func(f"{left}  {right}", f, [-.1, 1.1])

        print()
        print()
        p.show(file_name="piecewise_polynomial.html")


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
