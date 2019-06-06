# This file provides utilities for constructing polynomial interpolants.
# 
# The following objects are provided:
# 
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

# A generic Polynomial class that stores coefficients in monomial
# form. Provides numerically stable evaluation, derivative, and string
# operations for convenience.
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

    # Construct the polynomial that is the derivative of this polynomial.
    def derivative(self, d=1):
        if (d > 1): return self.derivative().derivative(d-1)
        else:       return Polynomial([c*i for (c,i) in zip(
                self.coefficients, range(len(self.coefficients)-1,0,-1))])

   # Given coefficients (and optionally points, assumed to be all zero)
    # for a Newton Polynomial, initialize a Monomial function.
    def __init__(self, coefficients, points=None):
        # If the user initialized this Polynomial with a Newton
        # Polynomial, then extract the points and coefficients.
        if (type(coefficients) == NewtonPolynomial):
            points = coefficients.points
            coefficients = coefficients.coefficients
        # If the user provided points, assume a Newton form polynomial
        # was provided and convert it to a monomial form.
        if (type(points) != type(None)):
            self.coefficients = to_monomial(coefficients, points)
        # Otherwise, assume standard polynomial coefficients were given.
        else:
            self.coefficients = coefficients

    # Evaluate this Polynomial at a point "x" in a numerically stable way.
    def __call__(self, x):
        total = self.coefficients[0]
        for d in range(1,len(self.coefficients)):
            total = self.coefficients[d] + x * total
        return total

    # Construct a string representation of this Polynomial.
    def __str__(self):
        s = ""
        for i in range(len(self.coefficients)):
            if   (i == len(self.coefficients)-1): x = ""
            elif (i == len(self.coefficients)-2): x = "x"
            else:   x = f"x^{len(self.coefficients)-1-i}"
            s += f"{self.coefficients[i]} {x}  +  "
        # Remove the trailing 
        return s[:-len("   +  ")]

# Extend the standard Polymomial class to hold Newton polynomials with
# points in addition to the coefficients.
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
# exactly match those values and derivatives at 0 on the left and 1
# on the right (any range can be mapped into this).
def polynomial_piece(left, right, interval=(0,1)):
    # Make sure both are lists.
    left  = list(left)
    right = list(right)
    # Rescale left and right to make their usage in the divided
    # difference table correct (by dividing by the factorial).
    mult = 1
    for i in range(2,len(left)):
        mult *= i
        left[i] /= mult
        right[i] /= mult
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
    return NewtonPolynomial(coefs, points)

# --------------------------------------------------------------------
#                            TESTING CODE

# Test the Polynomial class for basic operation.
def _test_Polynomial():
    f = Polynomial([3,2,1])
    assert(str(f) == "3 x^2  +  2 x  +  1")
    assert(str(f.derivative()) == "6 x  +  2")
    assert(str(f.derivative(2)) == "6")
    assert(str(f.derivative(3)) == "")
    assert(str(f.derivative(4)) == "")
    f = Polynomial([3,2,1],[0,0,0])
    assert(str(f) == "3 x^2  +  2 x  +  1")
    assert(str(f.derivative()) == "6 x  +  2")
    assert(str(f.derivative(2)) == "6")
    assert(str(f.derivative(3)) == "")
    assert(str(f.derivative(4)) == "")
    f = Polynomial([-1,10,-16,24,32,-32], [1,1,1,-1,-1,-1])
    assert(str(f) == "-1 x^5  +  9 x^4  +  6 x^3  +  -22 x^2  +  11 x  +  -3")
    assert(str(f.derivative()) == "-5 x^4  +  36 x^3  +  18 x^2  +  -44 x  +  11")

# Test the Polynomial class for basic operation.
def _test_NewtonPolynomial():
    f = NewtonPolynomial([-1,2], [1,-1])
    assert(str(f) == "2 + (x + 1)(-1)")
    assert(str(Polynomial(f)) == "-1 x  +  1")
    f = NewtonPolynomial([-1,10,-16,24,32,-32], [1,1,1,-1,-1,-1])
    assert(str(f) == "-32 + (x + 1)(32 + (x + 1)(24 + (x + 1)(-16 + (x - 1)(10 + (x - 1)(-1)))))")

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

if __name__ == "__main__":
    # Run the tests on this file.
    _test_Polynomial()
    _test_NewtonPolynomial()
    _test_polynomial()
    _test_polynomial_piece()
