# This file provides utilities for constructing polynomial interpolants.
# 
# The following objects are provided:
# 
#   Spline           -- A piecewise polynomial interpolant.
#   Polynomial       -- A monomial with stored coefficients,
#                       evaluation, derivative, and a string method.
#   NewtonPolynomial -- A Newton polynomial with stored coefficients,
#                       points, evaluation, derivative, and a string method.
# 
# The following functions are provided:
# 
#   polynomial -- Given x and y values, this produces a Newton
#                 polynomial that interpolates the provided points.
#   polynomial_piece -- Given a function value and derivatives,
#                       produce a polynomial that interpolates these
#                       values at the endpoints of an interval.
# 

# A piecewise polynomial function that supports evaluation,
# differentiation, and stitches together an arbitrary sequence of
# function values (and any number of derivatives) at data points
# (knots). Provides evaluation, derivation, and string methods.
# 
# Spline(knots, values):
#   Given a sequence of "knots" [float, ...] and an equal-length
#   sequence of "values" [[float, ...], ...] that define the
#   function value and any number of derivatives at every knot,
#   construct the piecewise polynomial that is this spline.
class Spline:
    # Define private internal variables for holding the knots, values,
    # and the functions. Provide access to "knots" and "values" properties.
    _knots = None
    _values = None
    _functions = None
    @property
    def knots(self): return self._knots
    @knots.setter
    def knots(self, knots): self._knots = list(knots)
    @property
    def values(self): return self._values
    @values.setter
    def values(self, values): self._values = [list(v) for v in values]
    
    def __init__(self, knots, values):
        assert(len(knots) == len(values))
        self.knots = knots
        self.values = values
        # Create the polynomial functions for all the pieces.
        self._functions = []
        for i in range(len(knots)-1):
            self._functions.append(
                polynomial_piece(self.values[i], self.values[i+1],
                                 (self.knots[i], self.knots[i+1]))
            )

    # Evaluate this Spline at a given x coordinate.
    def __call__(self, x):
        # Deterimine which interval this "x" values lives in.
        if (x <= self.knots[0]):    return self._functions[0](x)
        elif (x >= self.knots[-1]): return self._functions[-1](x)
        # Find the interval in which "x" exists.
        for i in range(len(self.knots)-1):
            if (self.knots[i] <= x <= self.knots[i+1]): break
        # If no interval was found, then something unexpected must have happened.
        else:
            class UnexpectedError(Exception): pass
            raise(UnexpectedError("This problem exhibited unexpected behavior."))
        # Now "self.knots[i] <= x <= self.knots[i+1]" must be true.
        return self._functions[i](x)
        
    # Compute the first derivative of this Spline.
    # WARNING: The returned spline does *not* have "values"!
    def derivative(self, d=1):
        # Create a spline (do not bother filling values) with the same
        # knot sequence and all the derivative functions.
        s = Spline([], [])
        s.knots = self.knots
        s._functions = [f.derivative(d) for f in self._functions]
        return s

    # Produce a string description of this spline.
    def __str__(self):
        s = "Spline:\n"
        s += f" [-inf, {self.knots[1]}]  =  "
        s += str(self._functions[0]) + "\n"
        for i in range(1,len(self.knots)-2):
            s += f" ({self.knots[i]}, {self.knots[i+1]}]  =  "
            s += str(self._functions[i]) + "\n"
        s += f" ({self.knots[-2]}, inf)  =  "
        s += str(self._functions[-1])
        return s

# A generic Polynomial class that stores coefficients in monomial
# form. Provides numerically stable evaluation, derivative, and string
# operations for convenience.
# 
# Polynomial(coefficients):
#    Given coefficients (or optionally a "NewtonPolynomial")
#    initialize this Monomial representation of a polynomial function.
class Polynomial:
    # Initialize internal storage for this Newton Polynomial.
    _coefficients = None
    # Protect the "coefficients" of this class with a getter and
    # setter to ensure a user does not break them on accident.
    @property
    def coefficients(self): return self._coefficients
    @coefficients.setter
    def coefficients(self, coefs): self._coefficients = list(coefs)
    # Define an alternative alias shortform "coefs".
    @property
    def coefs(self): return self._coefficients
    @coefs.setter
    def coefs(self, coefs): self._coefficients = list(coefs)

    def __init__(self, coefficients):
        # If the user initialized this Polynomial with a Newton
        # Polynomial, then extract the points and coefficients.
        if (type(coefficients) == NewtonPolynomial):
            coefficients = to_monomial(coefficients.coefficients,
                                       coefficients.points)
        self.coefficients = coefficients

    # Evaluate this Polynomial at a point "x" in a numerically stable way.
    def __call__(self, x):
        if (len(self.coefficients) == 0): return 0
        total = self.coefficients[0]
        for d in range(1,len(self.coefficients)):
            total = self.coefficients[d] + x * total
        return total

    # Construct the polynomial that is the derivative of this polynomial.
    def derivative(self, d=1):
        if (d == 0):  return self
        elif (d > 1): return self.derivative().derivative(d-1)
        else:         return Polynomial([c*i for (c,i) in zip(
                self.coefficients, range(len(self.coefficients)-1,0,-1))])

    # Construct a string representation of this Polynomial.
    def __str__(self):
        s = ""
        for i in range(len(self.coefficients)):
            if (self.coefficients[i] == 0): continue
            if   (i == len(self.coefficients)-1): x = ""
            elif (i == len(self.coefficients)-2): x = "x"
            else:   x = f"x^{len(self.coefficients)-1-i}"
            s += f"{self.coefficients[i]} {x}  +  "
        # Remove the trailing 
        return s.rstrip(" +")

# Extend the standard Polymomial class to hold Newton polynomials with
# points in addition to the coefficients.
# 
# NewtonPolynomial(coefficients, points):
#    Given a set of coefficients and a set of points (offsets), of
#    the same length, construct a standard Newton Polynomial.
class NewtonPolynomial(Polynomial):
    _points = None
    @property
    def points(self): return self._points
    @points.setter
    def points(self, points): self._points = list(points)

    # Store the coefficients and points for this Newton Polynomial.
    def __init__(self, coefficients, points):
        if (len(points) != len(coefficients)):
            raise(IndexError)
        self.coefficients = coefficients
        self.points = points

    # Construct the polynomial that is the derivative of this polynomial.
    def derivative(self, d=1): return Polynomial(self).derivative(d)

    # Evaluate this Newton Polynomial in a numerically stable way.
    def __call__(self, x):
        total = self.coefficients[0]
        for d in range(1,len(self.coefficients)):
            total = self.coefficients[d] + (x - self.points[d]) * total
        return total

    # Construct a string representation of this Newton Polynomial.
    def __str__(self):
        s = f"{self.coefficients[0]}"
        for i in range(1,len(self.coefficients)):
            sign = "-" if (self.points[i] >= 0) else "+"
            s = f"{self.coefficients[i]} + (x {sign} {abs(self.points[i])})({s})"
        return s

# Given Newton form coefficients and points, convert them to monomial
# form (where all points are 0) having only coefficients.
def to_monomial(coefficients, points):
    coefs = [coefficients[0]]
    for i in range(1,len(coefficients)):
        # Compute the old coefficients multiplied by a constant and
        # add the lower power coefficients that are shifted up.
        coefs.append(coefficients[i])
        coefs = [coefs[0]] + [coefs[j+1]-points[i]*coefs[j] for j in range(len(coefs)-1)]
    return coefs

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
    # Return the interpolating polynomial.
    return NewtonPolynomial(coefs, points)

# Given a (left value, left d1, ...), (right value, right d1, ...)
# pair of tuples, return the lowest order polynomial necessary to
# exactly match those values and derivatives at interval[0] on the left
# and interval[1] on the right ("interval" is optional, default [0,1]).
def polynomial_piece(left, right, interval=(0,1)):
    # Make sure both are lists.
    left  = list(left)
    right = list(right)
    # Fill values by matching them on both sides of interval (reducing order).
    for i in range(len(left) - len(right)):
        right.append( left[len(right)] )
    for i in range(len(right) - len(left)):
        left.append( right[len(left)] )
    # Rescale left and right to make their usage in the divided
    # difference table correct (by dividing by the factorial).
    mult = 1
    for i in range(2,len(left)):
        mult *= i
        left[i] /= mult
        right[i] /= mult
    # First match the function value, then compute the coefficients
    # for all of the higher order terms in the polynomial.
    coefs = list(left)
    # Compute the divided difference up until we flatten out the
    # highest provided derivative between the two sides.
    interval_width = interval[1] - interval[0]
    if (len(left) == 1): dds = []
    else:                dds = [(right[0] - left[0]) / interval_width]
    for i in range(len(left)-2):
        new_vals = [ (dds[0] - left[i+1]) / interval_width ]
        for j in range(len(dds)-1):
            new_vals.append( (dds[j+1] - dds[j]) / interval_width )
        new_vals.append( (right[i+1] - dds[-1]) / interval_width )
        dds = new_vals
    # Now the last row of the dd table should be level with "left" and
    # "right", we can complete the rest of the table in the normal way.
    row = [left[-1]] + dds + [right[-1]]
    # Build out the divided difference table.
    while (len(row) > 1):
        row = [(row[i+1]-row[i])/interval_width for i in range(len(row)-1)]
        coefs.append(row[0])
    # Reverse the coefficients to go from highest order (most nested)
    # to lowest order, making evaluation slightly more clean.
    points = [interval[1]]*len(left) + [interval[0]]*len(left)
    coefs = list(reversed(coefs))
    # Return the polynomial function.
    return Polynomial(to_monomial(coefs, points))

# --------------------------------------------------------------------
#                            TESTING CODE

# Test the Polynomial class for basic operation.
def _test_Polynomial():
    f = Polynomial([3,0,1])
    assert(str(f) == "3 x^2  +  1")
    f = Polynomial([3,2,1])
    assert(str(f) == "3 x^2  +  2 x  +  1")
    assert(str(f.derivative()) == "6 x  +  2")
    assert(str(f.derivative(2)) == "6")
    assert(str(f.derivative(3)) == "")
    assert(str(f.derivative(4)) == "")
    assert(f.derivative(3)(10) == 0)
    f = Polynomial(NewtonPolynomial([3,2,1],[0,0,0]))
    assert(str(f) == "3 x^2  +  2 x  +  1")
    assert(str(f.derivative()) == "6 x  +  2")
    assert(str(f.derivative(2)) == "6")
    assert(str(f.derivative(3)) == "")
    assert(str(f.derivative(4)) == "")
    assert(f.derivative(3)(5) == 0)
    f = Polynomial(to_monomial([-1,10,-16,24,32,-32], [1,1,1,-1,-1,-1]))
    assert(str(f) == "-1 x^5  +  9 x^4  +  6 x^3  +  -22 x^2  +  11 x  +  -3")
    assert(str(f.derivative()) == "-5 x^4  +  36 x^3  +  18 x^2  +  -44 x  +  11")

# Test the Polynomial class for basic operation.
def _test_NewtonPolynomial():
    f = NewtonPolynomial([-1,2], [1,-1])
    assert(str(f) == "2 + (x + 1)(-1)")
    assert(str(Polynomial(f)) == "-1 x  +  1")
    f = NewtonPolynomial([-1,10,-16,24,32,-32], [1,1,1,-1,-1,-1])
    assert(str(f) == "-32 + (x + 1)(32 + (x + 1)(24 + (x + 1)(-16 + (x - 1)(10 + (x - 1)(-1)))))")

# Test the "polynomial" interpolation routine (uses Newton form).
def _test_polynomial():
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

# Test the "polynomial_piece" interpolation routine.
def _test_polynomial_piece(plot=False):
    if plot:
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
        ([0,1,10], [0,-1,-10]),
        ([-2,2,10,6], [0,0,20,-6]),
    ]
    # Plot a bunch of sample functions.
    for (left, right) in left_rights:
        f = polynomial_piece( left, right )
        exact_coefs = list(map(round,f.coefs))[::-1]
        left_evals = [f(0)]
        right_evals = [f(1)]
        for i in range(1,len(left)):
            df = f.derivative(i)
            left_evals.append( df(0) )
            right_evals.append( df(1) )
        # TODO: Print out an error if the assert statement fails.
        assert(0 == sum(abs(true - app) for (true, app) in zip(left,left_evals)))
        assert(0 == sum(abs(true - app) for (true, app) in zip(right,right_evals)))
        # Create a plot of the functions if a demo is desired.
        if plot: p.add_func(f"{left}  {right}", f, [-.1, 1.1],
                            )#mode="markers", marker_size=2)
    if plot: p.show(file_name="piecewise_polynomial.html")

# Test the Spline class for basic operation.
def _test_Spline():
    knots = [0,1,2,3,4]
    values = [[0],[1,-1,0],[0,-1],[1,0,0],[0]]
    f = Spline(knots, values)
    for (k,v) in zip(knots,values):
        for d in range(len(v)):
            assert(f.derivative(d)(k) == v[d])

if __name__ == "__main__":
    # Run the tests on this file.
    _test_Polynomial()
    _test_NewtonPolynomial()
    _test_polynomial()
    _test_polynomial_piece()
    _test_Spline()
