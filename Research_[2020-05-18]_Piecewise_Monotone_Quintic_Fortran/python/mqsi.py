# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#                        fraction.py
# 
# Fraction class that provides infinite-precision rational numbers
# (making use of Python default arbitrarily large integers).

from math import gcd
class UsageError(Exception): pass

# This class-operator wrapper automatially returns a `float` type if
# the operand is neither the same type as the operator nor an `int`.
def float_fallback(operator):
    def wrapped_operator(primary, other):
        # Call the standard operator if the types match.
        if (type(other) == type(primary)):
            return operator(primary, other)
        # Handle integers in a special way (so that users can put ints
        # in code without breaking fractions).
        elif (type(other) == int):
            return operator(primary, type(primary)(other))
        # Otherwise, cast "other" to the type of "primary", call the
        # operator, then cast the result back to a 'float' type.
        else:
            return float(operator(primary, type(primary)(other)))
    return wrapped_operator

# This class-operator wrapper automatially casts to and returns an
# expected type (in common usage here, `primary` will be `Fraction`).
def auto_same_type(operator):
    def wrapped_operator(primary, other):
        # Call the standard operator if the types match.
        if (type(other) == type(primary)):
            return operator(primary, other)
        # Otherwise, cast "other" into the same type.
        else: return operator(primary, type(primary)(other))
    return wrapped_operator

# This can be 'float_fallback' or 'auto_same_type'.
TYPE_HANDLER = auto_same_type
ACCURACY = float('inf')

# A fraction instance. Initialized as:
# 
#     Fraction(numerator=0, denominator=1)
# 
# The arguments to the Fraction must either be Python integers or
# convertible to Python float via the 'float' built-in function.
class Fraction: 
    def __init__(self, numerator=0, denominator=1, reduce=True):
        # If "None" was provided, then set this fraction to have "None" values.
        if ((numerator is None) or (denominator is None)):
            self._numerator = None
            self._denominator = None
            return
        # If a Fraction was provided, shortcut and exit.
        elif (type(numerator) == Fraction):
            self._numerator   = numerator.numerator
            self._denominator = numerator.denominator
            return
        # Convert input arguments to the appropriate type.
        if (type(numerator) != int):
            try:    numerator = float(numerator)
            except: raise(UsageError(f"Unsupported numerator number type '{type(numerator)}', could not convert to float."))
        if (type(denominator) != int):
            try:    denominator = float(denominator)
            except: raise(UsageError(f"Unsupported denominator number type '{type(denominator)}', could not convert to float."))
        # Convert non-integer numerator and denominator appropriately.
        if (type(numerator) == float):
            # Handle 'infinity' in an expected way (without errors).
            if (abs(numerator) == float('inf')):
                if (abs(denominator) != float('inf')):
                    v = int(numerator < 0)
                    numerator, denominator = (-1)**(numerator < 0), 0
                elif (denominator < 0):
                    numerator, denominator = -numerator, -denominator
            # The normal usage case.
            else:
                numerator, _ = numerator.as_integer_ratio()
                denominator *= _
        if (type(denominator) == float):
            # Handle 'infinity' in an expected way (without errors).
            if (abs(denominator) == float('inf')):
                if (abs(numerator) != float('inf')):
                    numerator, denominator = 0, 1
            # The normal usage case.
            else:
                denominator, _ = denominator.as_integer_ratio()
                numerator *= _
        # Check to make sure the "infinity" cases are handled consistently.
        if (numerator == 0):
            if (denominator != 0): denominator //= abs(denominator)
        elif (denominator == 0):   numerator   //= abs(numerator)
        # Transfer the negativity to the numerator if possible.
        if   (numerator < 0) and (denominator < 0):
            numerator, denominator = -numerator, -denominator
        elif (numerator >= 0) and (denominator < 0):
            numerator, denominator = -numerator, -denominator
        # If this is indeterminate form, do not reduce!
        if (numerator == 0 == denominator): reduce = False
        # If the number has too much accuracy, round it.
        global ACCURACY
        if (denominator > ACCURACY):
            numerator = numerator * ACCURACY // denominator
            denominator = ACCURACY
        # Simplify the numerator and denominator if appropriate.
        if reduce:
            divisor = gcd(numerator, denominator)
            numerator   //= divisor
            denominator //= divisor
        # Store the numerator and denominator of this fraction.
        self._numerator = numerator
        self._denominator = denominator

    # Provide access to the numerator and denominator of this Fraction.
    @property
    def numerator(self): return self._numerator
    @property
    def denominator(self): return self._denominator

    # Convert this fraction into a float.
    def __float__(self):
        if (self.denominator == 0):
            if   (self.numerator > 0): return float('inf')
            elif (self.numerator < 0): return -float('inf')
            else:                      return float(0)
        else: return (self.numerator / self.denominator)
    # Convert this fraction to an integer.
    def __int__(self):
        if (0 == self.denominator == self.numerator): return 0
        return (self.numerator // self.denominator)
    # Provide a represenation of this Fraction.
    def __repr__(self): return f"{self.__class__.__name__}({self.numerator}, {self.denominator})"
    # Provide a string representation of this Fraction.
    def __str__(self):
        if (self.denominator == 0):
            if  (self.numerator == 0): return '0 / 0'
            elif (self.numerator > 0): return '+inf'
            elif (self.numerator < 0): return '-inf'
        elif (self.denominator == 1): return str(self.numerator)
        my_string = f"{self.numerator} / {self.denominator}"
        float_string = str(float(self))
        if (2*len(float_string) < len(my_string)): return float_string
        else:                                    return my_string
    # Format this Fraction by converting it to a float.
    def __format__(self, *args, **kwargs): return float(self).__format__(*args,**kwargs)
    # Create a "hash" of this object.
    def __hash__(self):
        if (float(self) == self): return hash(float(self))
        else:                     return hash(repr(self))
    # Add two fractions.
    @TYPE_HANDLER
    def __add__(a, b):
        if (0 == a.denominator == a.numerator): a = Fraction(0,1)
        if (0 == b.denominator == b.numerator): b = Fraction(0,1)
        an, bn = a.numerator, b.numerator
        ad, bd = a.denominator, b.denominator
        return Fraction( an * bd + bn * ad,  ad * bd )
    @TYPE_HANDLER
    def __radd__(b, a): return a + b

    #  Subtract two fractions.
    @TYPE_HANDLER
    def __sub__(a, b):
        if (0 == a.denominator == a.numerator): a = Fraction(0,1)
        if (0 == b.denominator == b.numerator): b = Fraction(0,1)
        an, bn = a.numerator, b.numerator
        ad, bd = a.denominator, b.denominator
        return Fraction( an * bd - bn * ad,  ad * bd )
    @TYPE_HANDLER
    def __rsub__(a, b): return -a + b

    # Multiply two fractions.
    @TYPE_HANDLER
    def __mul__(a, b):
        return Fraction(a.numerator * b.numerator, a.denominator * b.denominator)
    @TYPE_HANDLER
    def __rmul__(a, b): return a * b

    # Divide two fractions.
    @TYPE_HANDLER
    def __truediv__(a, b):
        return Fraction(a.numerator * b.denominator, a.denominator * b.numerator)
    @TYPE_HANDLER
    def __rtruediv__(b, a):
        return a / b

    # Divide two fractions (with integer result).
    def __floordiv__(a, b):
        b = Fraction(b)
        result = a / b
        result.numerator -= result.numerator % result.denominator
        return result.numerator // result.denominator
    def __rfloordiv__(b, a):
        a = Fraction(a)
        return a // b

    # Compute a mod b (and b mod a).
    @TYPE_HANDLER
    def __mod__(a, b):
        mult = a // b
        return a - b * mult
    @TYPE_HANDLER
    def __rmod__(b, a):
        return a % b

    # Compute a raised to the power b.
    # 
    # WARNING: Precision is lost if a root is taken (fractional power).
    def __pow__(a, b):
        # Handle the integer power case.
        if (type(b) == int):
            if (b == 0):  return Fraction(1,1,reduce=False)
            elif (b > 0): return Fraction(a.numerator**b, a.denominator**b)
            else:         return Fraction(a.denominator**(abs(b)), a.numerator**(abs(b)))
        else:
            # Handle the fractional power case (break up into two steps).
            b = Fraction(b)
            # Handle infinite values in a reasonable way.
            if (b.denominator == 0):
                # Raising to 0 power returns 1
                if (b.numerator == 0): return Fraction(1)
                # Raising to negative powers inverts the result.
                elif (b.numerator < 0): a = 1 / a
                # Handle balues based on "a" to determine where the limit is.
                if   (abs(a.numerator) > a.denominator): return Fraction(1,0)
                elif (abs(a.numerator) < a.denominator): return Fraction(0,1)
                else:                                    return a
            # If this is essentially an integer, just use that operation.
            elif (b.denominator == 1): return a ** b.numerator
            # First do the integer power (exact operation).
            intermediate = a ** abs(b.numerator)
            # Second do the inexact operation (converts to a float).
            result = float(intermediate) ** (1 / b.denominator)
            if (b.numerator < 0): result = 1 / result
            return result

    # When taking a Fraction as a power..
    def __rpow__(b, a):
        return Fraction(a) ** b

    # Negate this Fraction.
    def __neg__(a): return Fraction(-a.numerator, a.denominator, reduce=False)

    # Make this Fraction positive.
    def __abs__(a): return Fraction(abs(a.numerator), a.denominator, reduce=False)

    # a equals b. (special case is False for comparison with None)
    def __eq__(a, b):
        b = Fraction(b)
        return (a.numerator == b.numerator) and (a.denominator == b.denominator)

    # a less than b.
    def __lt__(a, b):
        b = Fraction(b)
        return a.numerator * b.denominator < a.denominator * b.numerator

    # a greater than b
    def __gt__(a, b):
        b = Fraction(b)
        return a.numerator * b.denominator > a.denominator * b.numerator

    # a less than or equal to b
    def __le__(a, b):
        b = Fraction(b)
        return a.numerator * b.denominator <= a.denominator * b.numerator

    # a greater than or equal to b
    def __ge__(a, b):
        b = Fraction(b)
        return a.numerator * b.denominator >= a.denominator * b.numerator

    # a not equal to zero
    def __bool__(a):
        return (a.numerator != 0) and (a.denominator != 0)

    # Return the rounded version of this Fraction.
    def __round__(self, *args, **kwargs):
        return round(float(self), *args, **kwargs)


# Some minor testing and demonstration code. Anything not tested here
# has not been explicitly tested!! Test cases will be built as more
# features are used.
def _test_Fraction(display=False, print=lambda *args, **kwargs: None):
    # If display is desired, re-assign the real print function.
    if display:
        import builtins
        print = builtins.print

    print()
    print("1/3              : ",1/3)
    print("a = Fraction(1,3): ",Fraction(1,3))
    print("b = Fraction(1/3): ",Fraction(1/3))

    a = Fraction(1,3)
    assert(a.numerator == 1)
    assert(a.denominator == 3)
    b = Fraction(1/3)
    assert(a != b)

    print()
    print("a - b:             ",a - b)
    assert((a - b) == Fraction(1, 54043195528445952))
    print("float(a - b):      ",float(a - b))
    print("a_sqrt = a**(1/2): ",a**(1/2))
    a_sqrt = a**(1/2)
    print()

    print("error = a - a_sqrt**2: ",a - a_sqrt**2)
    error = a - a_sqrt**2
    print("float(error):          ",float(error))
    print("error - float(a-b):    ",float(error - float(a-b)))
    # Error should certainly be less than SQRT error for rounded fractions.
    assert(error <= float(a-b)**(1/2))  
    print()

    print("Fraction():       ",Fraction())
    assert(Fraction() == Fraction(0,1))
    print("bool(Fraction()): ",bool(Fraction()))
    assert(not bool(Fraction()))
    print("bool(a):          ",bool(a))
    assert(bool(a))
    print()

    print("c = Fraction(float('inf')): ",Fraction(float('inf')))
    c = Fraction(float('inf'))
    assert(c == Fraction(1,0))
    print("Fraction(-float('inf')):    ",Fraction(-float('inf')))
    assert(Fraction(-float('inf')) == Fraction(-1,0))
    print("10 / c:    ",Fraction(10) / c)
    assert((Fraction(10) / c) == Fraction(0,1.))
    print("c / 10:    ",c / Fraction(10))
    assert((c / Fraction(10)) == Fraction(1.,0))
    print("c * 100.0: ",c * 100.0)
    assert((c * 100.0) == float('inf'))
    print("-c * 10:   ",-c)
    assert((-c*10) == Fraction(-1,0))
    print("c - c:     ",c - c)
    assert((c-c) == Fraction(0,0))
    print("float(c-c) ", float(c-c))
    assert(float(c-c) == 0.0)
    print()
    print("a * (c - c)  ", a*(c - c))
    assert((a*(c-c)) == Fraction(0,0))
    print("a / (c - c)  ", a/(c - c))
    assert((a/(c-c)) == Fraction(0,0))
    print("a ** c       ", a**c)
    assert((a**c) == Fraction(0,1))
    print("a ** (-c)    ", a**(-c))
    assert((a**(-c)) == Fraction(1,0))
    print("(-a) ** c    ", (-a)**c)
    assert(((-a)**c) == Fraction(0,1))
    print("(-a) ** (-c) ", (-a)**(-c))
    assert(((-a)**(-c)) == Fraction(1,0))
    print()

    print("d = Fraction(1,float('inf')): ",Fraction(1,float('inf')))
    d = Fraction(1,float('inf'))
    assert(d == Fraction(0,1))
    print("Fraction(1,-float('inf')):    ",Fraction(1,-float('inf')))
    assert(Fraction(1,-float('inf')) == Fraction(0,1))
    print("10 / d:  ",10 / d)
    assert((10/d) == Fraction(1,0))
    print("d / 10:  ",d / 10)
    assert((d/10) == Fraction(0,1))
    print("d * 100: ",d * 100)
    assert((d * 100) == Fraction(0,1))
    print("c * d:   ",c * d)
    assert((c * d) == Fraction(0,0))
    print("-d:      ",-d)
    assert((-d) == Fraction(0,1))
    print()

    # Make sure the "Fraction" object is uniquely hashable.
    assert(hash(a) == hash(Fraction(1,3)))
    assert(hash(b) == hash(Fraction(1/3)))
    assert(hash(c) == hash(1/d))

    # Make sure that the power operator behaves in expected ways.
    a = Fraction(7,3)
    assert(a**2 == a**2.)
    assert(a**(-2) == a**(-2.))
    assert(Fraction(2.0, -1/2) == Fraction(-4))

# Convert a value or list of values into a list of Fraction objects.
def make_exact(value):
    try:    return list(map(make_exact, value))
    except: return Fraction(value)

if __name__ == "__main__":
    print()
    print("Testing fraction.py")
    _test_Fraction(display=False)
    print(" all PASSED.")


# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#                        polynomial.py
# 
# This file provides utilities for constructing polynomial interpolants.
# 
# The following objects are provided:
# 
#   Spline           -- A piecewise polynomial interpolant with
#                       evaluation, derivative, integration, negation,
#                       multiplication, addition, and a string method.
#   Polynomial       -- A polynomial with nonrepeating coefficients,
#                       evaluation, derivative, and a string method.
#   NewtonPolynomial -- A Newton polynomial with stored coefficients,
#                       offsets, evaluation, derivative, and a string method.
# 
# The following functions are provided:
# 
#   fit              -- Use a local polynomial interpolant to estimate
#                       derivatives and construct an interpolating
#                       Spline with a specified level of continuity.
#   polynomial       -- Given x and y values, this produces a minimum
#                       degree Newton polynomial that interpolates the
#                       provided points (this is a *single* polynomial).
# 

# This general purpose exception will be raised during user errors.
class UsageError(Exception): pass

# A piecewise polynomial function that supports evaluation,
# differentiation, and stitches together an arbitrary sequence of
# function values (and any number of derivatives) at data points
# (knots). Provides evaluation, derivation, and string methods.
# It is recommended to use EXACT ARITHMETIC internally (which will
# automatically return floats unless Fraction objects are provided as
# evaluation points). Exact arithmetic is achieved by providing knots
# and values composed of Fraction objects.
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
    _derivative = 0
    @property
    def knots(self): return self._knots
    @knots.setter
    def knots(self, knots): self._knots = make_exact(knots)
    @property
    def values(self): return self._values
    @values.setter
    def values(self, values): self._values = make_exact(values)
    @property
    def functions(self): return self._functions
    @functions.setter
    def functions(self, functions): self._functions = list(functions)
    
    # Create a piecewise polynomial that matches the given values at
    # the given knots, alternatively provide a list of functions that
    # exist over a set of intervals defined by knots.
    def __init__(self, knots, values=None, functions=None):
        assert(len(knots) >= 1)
        self.knots = knots
        # Store the 'values' at each knot if they were provided.
        if (values is not None):
            assert(len(knots) == len(values))
            self.values = values
        # Store the 'functions' over each interval if they were provided.
        if (functions is not None):
            # TODO: Verify that the provided functions match the values.
            assert(len(functions) == len(knots)-1)
            self.functions = functions
        # Use the 'values' to generate 'functions' if no functions were provided.
        elif (values is not None):
            # Create the polynomial functions for all the pieces.
            self._functions = []
            for i in range(len(knots)-1):
                v0, v1 = self.values[i], self.values[i+1]
                k0, k1 = self.knots[i], self.knots[i+1]
                x = [k0, k1]
                y = [v0[0], v1[0]]
                derivs = {i*'d'+'x0':v0[i] for i in range(1,len(v0))}
                derivs.update({i*'d'+'x1':v1[i] for i in range(1,len(v1))})
                # Construct a local polynomial that interpolates the 
                # function (and derivative) values.
                f = polynomial(x, y, **derivs)
                # Store the function (assuming it's correct).
                self._functions.append( f )
        # No 'values' nor 'functions' were provided, usage error.
        else: raise(UsageError("Either 'values' or 'functions' must be provided with knots to construct a spline."))

    # Evaluate this Spline at a given x coordinate.
    def function_at(self, x):
        # If "x" was given as an iterable, then iterate over it.
        try:    return [self.function_at(v) for v in x]
        except: pass
        # Deterimine which interval this "x" values lives in.
        if (x <= self.knots[0]):    return self._functions[0]
        elif (x >= self.knots[-1]): return self._functions[-1]
        # Find the interval in which "x" exists.
        for i in range(len(self.knots)-1):
            if (self.knots[i] <= x < self.knots[i+1]): break
        # If no interval was found, then something unexpected must have happened.
        else:
            class UnexpectedError(Exception): pass
            raise(UnexpectedError("This problem exhibited unexpected behavior. Check code."))
        # Return the applicable function.
        return self._functions[i]

    # Compute the integral of this Spline.
    def integral(self, i=1, c=None): return self.derivative(-i, c)

    # Compute the first derivative of this Spline.
    def derivative(self, d=1, c=None):
        # For integration, adjust the additive term to reflect the
        # expected lefthand side value of each function.
        if (d < 0):
            deriv_funcs = self._functions
            for i in range(-d):
                deriv_funcs = [f.integral(1) for f in deriv_funcs]
                if (c is None): total = self(self.knots[0])
                else:           total = c
                for i in range(len(deriv_funcs)):
                    deriv_funcs[i].coefficients[-1] = (
                        total - deriv_funcs[i](self.knots[i]))
                    total = deriv_funcs[i](self.knots[i+1])
        else:
            # Create a spline with the same knot sequence and all the
            # derivative functions and associated values.
            deriv_funcs = self._functions
            for i in range(d):
                deriv_funcs = [f.derivative(1) for f in deriv_funcs]

        # Construct the new spline, pass the "values" even though
        # nothing will be done with them. Assign the new "functions".
        s = Spline(self.knots, self.values, functions=deriv_funcs)
        s._derivative = self._derivative + d
        # Return the new derivative Spline.
        return s

    # Evaluate this Spline at a given x coordinate.
    def __call__(self, x):
        # If "x" was given as an iterable, then iterate over it.
        try:    return [self(v) for v in x]
        except: pass
        # Get the appropriate function and compute the output value.
        return self.function_at(x)(x)
        
    # Add this "Spline" object to another "Spline" object.
    def __add__(self, other):
        # Check for correct usage.
        if (type(other) != type(self)):
            raise(UsageError(f"Only '{type(self)} objects can be added to '{type(self)}' objects, but '{type(other)}' was given."))
        # Generate the new set of knots.
        knots = sorted(set(self._knots + other._knots))
        # Compute the functions over each interval.
        functions = []
        for i in range(len(knots)-1):
            # Get the knot, nearby knots, and order of resulting
            # polynomial at this particular knot.
            left, right = knots[i], knots[i+1]
            k = knots[i]
            my_poly = self.function_at(k)
            other_poly = other.function_at(k)
            order = max(len(my_poly.coefficients), len(other_poly.coefficients))
            # Evaluate the function at equally spaced "x" values, TODO:
            # this should be Chebyshev nodes for numerical stability.
            x = [(step / (order-1)) * (right - left) + left for step in range(order)]
            y = [self(node) + other(node) for node in x]
            # Construct the interpolating polynomial.
            functions.append( polynomial(x, y) )
            f = functions[-1]
        # Return the new added function.
        return Spline(knots, functions=functions)

    # Multiply this "Spline" obect by another "Spline".
    def __mul__(self, other):
        # Check for correct usage.
        if (type(other) != type(self)):
            raise(UsageError(f"Only '{type(self)} objects can be multiplied by '{type(self)}' objects, but '{type(other)}' was given."))
        # Generate the new set of knots.
        knots = sorted(set(self._knots + other._knots))
        # Compute the functions over each interval.
        functions = []
        for i in range(len(knots)-1):
            # Get the knot, nearby knots, and order of resulting
            # polynomial at this particular knot.
            left, right = knots[i], knots[i+1]
            k = knots[i]
            my_poly = self.function_at(k)
            other_poly = other.function_at(k)
            order = max(2,len(my_poly.coefficients) + len(other_poly.coefficients) - 1)
            # Evaluate the function at equally spaced "x" values, TODO:
            # this should be Chebyshev nodes for numerical stability.
            x = [(step / (order-1)) * (right - left) + left for step in range(order)]
            y = [self(node) * other(node) for node in x]
            # Construct the interpolating polynomial.
            functions.append( polynomial(x, y) )
            f = functions[-1]
        # Return the new added function.
        return Spline(knots, functions=functions)

    # Raise an existing spline to a power.
    def __pow__(self, number):
        if (type(number) != int) or (number < 1):
            raise(TypeError(f"Only possible to raise '{type(self)}' to integer powers greater than  or equal to 1."))
        # Start with a copy of "self", multiply in until complete.
        outcome = Spline(self.knots, values=self.values,
                         functions=[Polynomial(f) for f in self.functions])
        for i in range(number-1): outcome = outcome * self
        return outcome

    # Subtract another spline from this spline (add after negation).
    def __sub__(self, other): return self + (-other)

    # Negate this spline, create a new one that is it's negative.
    def __neg__(self):
        # Create a new spline, but negate all internal function coefficients.
        return Spline(self.knots, functions=[-f for f in self.functions])

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

# A generic Polynomial class that stores coefficients in a fully
# expanded form. Provides numerically stable evaluation, derivative,
# integration, and string operations for convenience. Coefficients
# should go highest to lowest order.
# 
# Polynomial(coefficients):
#    Given coefficients (or optionally a "NewtonPolynomial")
#    initialize this Polynomial representation of a polynomial function.
# 
# EXAMPLE:
# Polynomial([a,b,c])
#   = a x^2  +  b x  +  c
class Polynomial:
    # Initialize internal storage for this Polynomial.
    _coefficients = None
    # Protect the "coefficients" of this class with a getter and
    # setter to ensure a user does not break them on accident.
    @property
    def coefficients(self): return self._coefficients
    @coefficients.setter
    def coefficients(self, coefs): self._coefficients = make_exact(coefs)
    # Initialize this Polynomial given coefficients.
    def __init__(self, coefficients=(0,)):
        # If the user initialized this Polynomial with a Newton
        # Polynomial, then extract the points and coefficients.
        if (type(coefficients) == NewtonPolynomial):
            # Unpack the NewtonPolynomial into coefficients and points
            # by multiplying everything out (and unrolling it).
            newton_polynomial = coefficients
            c, p = newton_polynomial.coefficients, newton_polynomial.zeros
            coefficients = [c[0]]
            for i in range(1,len(c)):
                # Compute the old coefficients multiplied by a constant and
                # add the lower power coefficients that are shifted up.
                coefficients.append(c[i])
                coefficients = [coefficients[0]] + \
                               [coefficients[j+1]-p[i]*coefficients[j]
                                for j in range(len(coefficients)-1)]
        # If the user initialized this Polynomial, do a copy.
        elif (type(coefficients) == type(self)):
            coefficients = coefficients.coefficients
        # Store the coeficients.
        self.coefficients = coefficients # <- implicit copy
        # Remove all leading 0 coefficients.
        for i in range(len(self._coefficients)):
            if (self._coefficients[i] != 0): break
        else: i = len(self._coefficients)-1 # <- include the last 0, if needed
        self._coefficients = self._coefficients[i:]

    # Evaluate this Polynomial at a point "x" in a numerically stable way.
    def __call__(self, x):
        try: return [self(v) for v in x]
        except TypeError: pass # <- don't worry if it is not iterable
        if (len(self._coefficients) == 0): return 0
        total = self._coefficients[0]
        for d in range(1,len(self._coefficients)):
            total = self._coefficients[d] + x * total
        return total

    # Construct the polynomial that is the integral of this polynomial.
    def integral(self, i=1): return self.derivative(-i)

    # Construct the polynomial that is the derivative (or integral) of this polynomial.
    def derivative(self, d=1):
        if (d == 0): return self
        elif (d > 1):  return self.derivative(1).derivative(d-1)
        elif (d == 1): return Polynomial([c*i for (c,i) in zip(
                self._coefficients, range(len(self._coefficients)-1,0,-1))])
        elif (d < -1):  return self.derivative(-1).derivative(d+1)
        elif (d == -1): return Polynomial([c/(i+1) for (c,i) in zip(
                self._coefficients, range(len(self._coefficients)-1,-1,-1))]+[0])

    # Determines if two polynomials are equivalent.
    def __eq__(self, other):
        if type(other) == NewtonPolynomial: other = Polynomial(other)
        return all(c1 == c2 for (c1, c2) in zip(self._coefficients, other._coefficients))

    # Negate this polynomial.
    def __neg__(self): return Polynomial([-c for c in self._coefficients])

    # Construct a string representation of this Polynomial.
    def __str__(self):
        s = ""
        for i in range(len(self._coefficients)):
            if (self._coefficients[i] == 0): continue
            if   (i == len(self._coefficients)-1): x = ""
            elif (i == len(self._coefficients)-2): x = "x"
            else:   x = f"x^{len(self._coefficients)-1-i}"
            s += f"{str(self._coefficients[i])} {x}  +  "
        # Remove the trailing 
        s = s.rstrip(" +")
        # Return the final string.
        return s

# Extend the standard Polymomial class to hold Newton polynomials with
# zeros in addition to the coefficients. This is more convenient when
# constructing interpolating polynomials from divided difference tables.
# 
# NewtonPolynomial(coefficients, zeros):
#    Given a set of coefficients and a set of zeros (offsets), of
#    the same length, construct a standard Newton
#    Polynomial. Coefficients are stored from highest order term to
#    lowest order. Earlier zeros are evaluated earlier.
# 
# EXAMPLE:
# NewtonPolynomial([a,b,c], [s1, s2, s3])
#   =  c + (x - s3)(b + (x - s2)(a))
# 
# NOTE:
#   Notice that "s1" is never used in the computation, as it is
#   redundant with the term "a".
class NewtonPolynomial(Polynomial):
    _zeros = None
    @property
    def zeros(self): return self._zeros
    @zeros.setter
    def zeros(self, zeros): self._zeros = make_exact(zeros)

    # Store the coefficients and zeros for this Newton Polynomial.
    def __init__(self, coefficients, zeros):
        self.coefficients = coefficients
        # Strip off the 0-valued coefficients (all these will be 0'd out).
        for i in range(len(self._coefficients)):
            if (self._coefficients[i] != 0): break
        else: i = len(self._coefficients)-1
        # Store locally.
        self._coefficients = self._coefficients[i:]
        self.zeros = [p for (j,p) in enumerate(zeros) if j >= i]

    # Construct the polynomial that is the derivative of this
    # polynomial by converting to polynomial form and differntiating.
    def derivative(self, d=1): return Polynomial(self).derivative(d)

    # Evaluate this Newton Polynomial (in a numerically stable way).
    def __call__(self, x):
        total = self._coefficients[0]
        for d in range(1,len(self._coefficients)):
            total = self._coefficients[d] + (x - self._zeros[d]) * total
        return total

    # Negate this polynomial.
    def __neg__(self): return -Polynomial(self)

    # Construct a string representation of this Newton Polynomial.
    def __str__(self):
        s = f"{str(self._coefficients[0])}"
        for i in range(1,len(self._coefficients)):
            sign = "-" if (self._zeros[i] >= 0) else "+"
            s = f"{str(self._coefficients[i])} + (x {sign} {str(abs(self._zeros[i]))})({s})"
        return s


# Given data points "x" and data values "y", construct an
# interpolating spline over the given points with specified level of
# continuity using a sufficiently continuous polynomial fit over
# neighboring points.
#  
# x: A strictly increasing sequences of numbers.
# y: Function values associated with each point.
# 
# continuity:
#   The level of continuity desired in the interpolating function.
# 
def fit(x, y, continuity=0):
    assert( len(x) == len(y) )
    # Sort the "x" values if they were not given in sorted order.
    if not all(x[i] < x[i+1] for i in range(len(x)-1)):
        indices = sorted(range(len(x)), key=lambda i: x[i])
        x = [x[i] for i in indices]
        y = [y[i] for i in indices]
    # Get the knots and initial values for the spline.
    knots = [v for v in x]
    values = [[v] for v in y]
    # Compute the local order of the spline that will be used to
    # interpolate the data (this is more than necessary, but makes
    # the final result much smoother than the lowest possible order).
    # In general, the highest order necessary is "continuity + 2".
    order = 2*(continuity+1)
    step = (order+1)//2
    # Construct local polynomial inteprolants to estimate the derivatives.
    for i in range(len(x)):
        lower = max(0, i-step)
        upper = min(len(x)-1, i+step)
        df = polynomial(x[lower:upper+1], y[lower:upper+1])
        for d in range(continuity):
            df = df.derivative()
            values[i].append(df(x[i]))
    # Return the interpolating spline.
    return Spline(knots, values)


# Given unique "x" values and associated "y" values (of same length),
# construct an interpolating polynomial with the Newton divided
# difference method. Allows for the incorporation of derivative
# constraints via keyword arguments. I.e. the keyword argument "dx0=0"
# makes the derivative of the polynomial 0 at the first element of x.
# Returns a NewtonPolynomial object.
# 
# EXAMPLE:
# polynomial([1, 2, 3], [-1, 2, -3], dx1=0)
#   => cubic polynomial interpolating the points (0,-1), (1,2), and
#      (2,-3) that has a first derivative of 0 at x=2.
def polynomial(given_x, given_y, **derivs):
    # Sort the data by "x" value.
    indices = sorted(range(len(given_x)), key=lambda i: given_x[i])
    # Construct the initial x and y (with repetitions to match derivatives).
    x = []; y = []; index_ranges = {}
    for i in indices:
        dxs = sorted(key for key in derivs if (i == int(key.split('x')[1])))
        largest_d = ([key.count("d") for key in dxs] + [0])[0]
        index_ranges[(len(x), len(x)+largest_d)] = i
        x += [given_x[i]] * (largest_d+1)
        y += [given_y[i]] * (largest_d+1)
    # Compute the divided difference table.
    divisor = 1
    dd_values = [y]
    for d in range(1, len(x)):
        divisor *= d
        slopes = []
        for i in range(len(dd_values[-1])-1):
            # Identify which "original point index" this spot belongs to.
            idx = ''
            for (start, end) in index_ranges:
                if (start <= i <= end-d): idx = index_ranges[(start,end)]
            # Substitute in derivative value if it was provided.
            key = d*"d"+f"x{idx}"
            if key in derivs: dd = derivs[key] / divisor
            # Otherwise a value wasn't provided, compute the divided difference.
            elif (x[i+d] == x[i]): dd = 0
            else: dd = (dd_values[-1][i+1] - dd_values[-1][i]) / (x[i+d] - x[i])
            slopes.append( dd )
        # Add in the finished next row of the divided difference table.
        dd_values.append( slopes )
    # Get the divided difference (polynomial coefficients) in reverse
    # order so that the most nested value (highest order) is first.
    # Return as an interpolating polynomial in Newton form.
    return NewtonPolynomial(
        (row[0] for row in reversed(dd_values)),  reversed(x))


# --------------------------------------------------------------------
#                            TESTING CODE

# Test the Spline class for basic operation.
def _test_Spline():
    knots = [0,1,2,3,4]
    values = [[0],[1,-1,0],[0,-1],[1,0,0],[0]]
    # Create the knots and values.
    knots = [Fraction(k) for k in knots]
    values = [[Fraction(v) for v in vals] for vals in values]
    # Create the spline.
    f = Spline(knots, values)
    for (k,v) in zip(knots,values):
        for d in range(len(v)):
            try: assert(f.derivative(d)(k) == v[d])
            except:
                print()
                print('-'*70)
                print("      TEST CASE")
                print("Knot:           ", k)
                print("Derivative:     ", d)
                print("Expected value: ", v[d])
                print("Received value: ", f.derivative(d)(k))
                print()
                print(f)
                print('-'*70)
                raise(Exception("Failed test case."))
    # TODO: Test cases for the addition, multiplication, negation
    #       and integration of a Spline object.
    # 
    # from util.plot import Plot
    # g = Spline(knots, values)
    # s_add = f + g
    # add = lambda x: f(x) + g(x)
    # s_mult = f * g
    # mult = lambda x: f(x) * g(x)
    # p = Plot()
    # p.add_func("f", f, [min(knots), max(knots)])
    # p.add_func("f'", f.derivative(1), [min(knots), max(knots)])
    # p.add_func("f''", f.derivative(2), [min(knots), max(knots)])
    # p.add_func("g", g, [min(knots), max(knots)])
    # p.add_func("s add", s_add, [min(knots), max(knots)])
    # p.add_func("add", add, [min(knots), max(knots)])
    # p.add_func("s mult", s_mult, [min(knots), max(knots)])
    # p.add_func("mult", mult, [min(knots), max(knots)])
    # # Negated
    # p.add_func("-f", -f, [min(knots), max(knots)])
    # p.add_func("-f true", lambda x: -(f(x)), [min(knots), max(knots)])    
    # p.add_func("-mult", -s_mult, [min(knots), max(knots)])
    # p.add_func("-mult true", lambda x: -s_mult(x), [min(knots), max(knots)])
    # p.add_func("int(f)", f.derivative(-1), [min(knots), max(knots)])
    # p.add_func("int(-mult)", (-s_mult).derivative(-1), [min(knots), max(knots)])
    # p.show()
    # exit()


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
    f = Polynomial(NewtonPolynomial([-1,10,-16,24,32,-32], [1,1,1,-1,-1,-1]))
    assert(str(f) == "-1 x^5  +  9 x^4  +  6 x^3  +  -22 x^2  +  11 x  +  -3")
    assert(str(f.derivative()) == "-5 x^4  +  36 x^3  +  18 x^2  +  -44 x  +  11")
    # Check that integrals work too.
    assert(str(f.derivative().derivative(-1)) == "-1 x^5  +  9 x^4  +  6 x^3  +  -22 x^2  +  11 x")
    assert(str(f.derivative().derivative(-1).derivative()) == "-5 x^4  +  36 x^3  +  18 x^2  +  -44 x  +  11")


# Test the Polynomial class for basic operation.
def _test_NewtonPolynomial():
    f = NewtonPolynomial([-1,2], [1,-1])
    assert(str(f) == "2 + (x + 1)(-1)")
    assert(str(Polynomial(f)) == "-1 x  +  1")
    f = NewtonPolynomial([-1,10,-16,24,32,-32], [1,1,1,-1,-1,-1])
    assert(str(f) == "-32 + (x + 1)(32 + (x + 1)(24 + (x + 1)(-16 + (x - 1)(10 + (x - 1)(-1)))))")


# Test the "polynomial" interpolation routine (uses Newton form).
def _test_polynomial(plot=True):
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


# Test the "fit" function. (there is testing code built in, so this
# test is strictly for generating a visual to verify).
def _test_fit(plot=False):
    x_vals = list(map(Fraction, [0,.5,2,3.5,4,5.3,6]))
    y_vals = list(map(Fraction, [1,2,2.2,3,3.5,4,4]))
    f = fit(x_vals, y_vals, continuity=2)
    # Execute with different operational modes, (tests happen internally).
    f = fit(x_vals, y_vals, continuity=2)
    for i in range(len(f._functions)):
        f._functions[i] = Polynomial(f._functions[i])
    if plot:
        print("f: ",f)
        from util.plot import Plot
        plot_range = [min(x_vals)-.1, max(x_vals)+.1]
        p = Plot()
        p.add("Points", list(map(float,x_vals)), list(map(float,y_vals)))
        p.add_func("f)", f, plot_range)
        p.add_func("f'", f.derivative(1), plot_range, dash="dash")
        p.add_func("f'')", f.derivative(2), plot_range, dash="dot")
        # def L2(f, low, upp):
        #     f2 = f**2
        #     p.add_func("f''<sup>2</sup>", f2, plot_range)
        #     int_f2 = f2.integral()
        #     return int_f2(upp) - int_f2(low)
        # print("L2(f''):", float(L2(f.derivative(2), x_vals[0], x_vals[-1])))
        p.show()


if __name__ == "__main__":
    # Run the tests on this file.
    print()
    print("Testing polynomial.py")
    print(" Spline")
    _test_Spline()
    print(" Polynomial")
    _test_Polynomial()
    print(" NewtonPolynomial")
    _test_NewtonPolynomial()
    print(" polynomial")
    _test_polynomial()
    print(" fit")
    _test_fit()
    print(" all PASSED.")


# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#                        monotone.py
# 
# Codes for creating a monotone quintic interpolating spline.
# 

# Given "i" the index at which the first derivative should be
# estimated, construct local quadratic fits and pick the slope of the
# one with the lowest curvature.
def quadratic_facet(x, y, i, extremes, flats):
    # If this is a local flat, estimate 0.
    if (i in flats): return (0, 0)
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
    if (len(derivatives) == 0): return (0, 0)
    # Sort the derivatives by magnitude of curvature, then return the 
    # values from the function with the least curvature.
    derivatives.sort(key=lambda d: abs(d[1]))
    return derivatives[0]


# Where "d" is the number of derivatives to estimate.
def weighted_harmonic(x, y, i, extremes, flats, d=2):
    if (d == 1):
        if (i == 0) or (i == len(x)-1):
            d_value = quadratic_facet(x,y,i,extremes,flats)[0]
        else:
            h1 = (x[i] - x[i-1])
            h2 = (x[i+1] - x[i])
            d1 = (y[i] - y[i-1]) / h1
            d2 = (y[i+1] - y[i]) / h2
            w1 = h1 + 2*h2
            w2 = 2*h1 + h2
            if (d1*d2 < 0):    d_value = 0
            elif (d1*d2 == 0): d_value = 0
            else:              d_value = (w1 + w2) / (w1/d1 + w2/d2)
        return (d_value,)
    else:
        if (i == 0) or (i == len(x)-1):
            d_values = quadratic_facet(x,y,i,extremes,flats) + (0,)*(max(0,d-2))
        else:
            d_left = weighted_harmonic(x, y, i-1, extremes, flats, d=d-1)
            d_center = weighted_harmonic(x, y, i, extremes, flats, d=d-1)
            d_right = weighted_harmonic(x, y, i+1, extremes, flats, d=d-1)
            h1 = (x[i] - x[i-1])
            h2 = (x[i+1] - x[i])
            d1 = (d_center[-1] - d_left[-1]) / h1
            d2 = (d_right[-1] - d_center[-1]) / h2
            w1 = h1 + 2*h2
            w2 = 2*h1 + h2
            if (d1*d2 < 0):    d_value = 0
            elif (d1*d2 == 0): d_value = 0
            else:              d_value = (w1 + w2) / (w1/d1 + w2/d2)
            d_values = d_center + (d_value,)
        return d_values


# Compute a monotone quintic spline that interpolates given data.
def monotone_quintic_spline(x, y, accuracy=2**(-26), verbose=False,
                            monotone=True, estimator=quadratic_facet):
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
    for i in range(len(x)): values[i] += estimator(
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

# Function for computing tight condition on monotonicity. Returns True
# if an interval is monotone, False otherwise.
def is_monotone(U0, U1, F0, DF0, DDF0, F1, DF1, DDF1, EPS=2**(-26)):
    # Flip the sign to only consider the monotone increasing case.
    if (F1 < F0):
        F0,DF0,DDF0 = -F0, -DF0, -DDF0
        F1,DF1,DDF1 = -F1, -DF1, -DDF1
    # Check for a "flat", if it is flat then only allow zeros.
    if (abs(F1-F0) < EPS*(1 + abs(F0) + abs(F1))):
        return (DF0 == 0) and (DF1 == 0) and (DDF0 == 0) and (DDF1 == 0)
    # Make sure the slopes point in the right direction.
    if (((F1 - F0) * DF0) < 0): return False
    if (((F1 - F0) * DF1) < 0): return False
    W = U1 - U0
    # Simplified cubic monotone case.
    if ((abs(DF0) < EPS) or (abs(DF1) < EPS)):
        if (DDF1*W > 4*DF1): return False
        TEMP = DF0 * (4*DF1 - DDF1*W)
        if (TEMP > 0): TEMP = 2 * TEMP**(0.5)
        if (TEMP + 3*DF0 + DDF0*W < 0): return False
        if (60*(F1-F0) - W*((24*DF0 + 32*DF1) - 2*TEMP +
                            W*(3*DDF0 - 5*DDF1)) < 0): return False
        return True
    # Full quintic monotone case.
    if (W*(2*(DF0*DF1)**(0.5) - 3*(DF0+DF1)) + 24*(F1-F0)) <= 0: return False
    TEMP = (DF0*DF1)**(0.75)
    ALPHA = (4*DF1 - DDF1*W) * (DF0)**(0.5) / TEMP
    GAMMA = (4*DF0 + DDF0*W) * (DF1)**(0.5) / TEMP
    BETA = (3*((DDF1-DDF0)*W - 8*(DF0+DF1)) + 60*(F1-F0)/W) / (2 * (DF0*DF1)**(0.5))
    if (BETA <= 6): TEMP = -(BETA + 2) / 2
    else:           TEMP = -2*(BETA - 2)**(0.5)
    return (ALPHA > TEMP) and (GAMMA > TEMP)

# Construct an interpolant that minimizes the L2 of the second derivative.
def min_curve(x, y, t=10):
    # Pick first and second derivative values that minimize the L2.
    import numpy as np
    from util.optimize import minimize, zero_on_line
    # Define an objective function that minimizes L2 and makes the derivative positive.
    def obj_func(v):
        f = Spline(x, list(zip(y,v[0::2],v[1::2])))
        df = f.derivative()
        dfx = df(x)
        ddf = df.derivative()
        ddfx = ddf(x)
        if not all(is_monotone(x[i],x[i+1],y[i],dfx[i],ddfx[i],
                               y[i+1],dfx[i+1],ddfx[i+1])
                   for i in range(len(x)-1)):
            return float('inf')
        L2 = ((ddf)**2).integral()
        return L2(x[-1]) - L2(x[0])
    # Get an initial solution by using the monotone quintic.
    d_vec = []
    f = monotone_quintic_spline(x,y,estimator=weighted_harmonic)
    df = f.derivative()
    ddf = df.derivative()
    for v in x: d_vec += [df(v), ddf(v)]
    # Try and get an existing solution.
    try:    from optimization_checkpoint import solution
    except: solution = None
    # Use the existing solution if it looks right.
    if (solution is not None) and (len(solution) == len(d_vec)):
        d_vec = solution
    # Get the best interpolating spline by minimize L2.
    v = minimize(obj_func, d_vec, max_time=t)
    return Spline(x, list(zip(y,v[0::2],v[1::2])))

# nonmonotone_funcs = []
# def is_monotone(U0, U1, F0, DF0, DDF0, F1, DF1, DDF1, **kwargs):
#     f = polynomial([U0,U1], [F0,F1], dx0=DF0, ddx0=DDF0, dx1=DF1, ddx1=DDF1)
#     df = f.derivative()
#     output = _is_monotone(U0, U1, F0, DF0, DDF0, F1, DF1, DDF1, **kwargs)
#     if not output:
#         print()
#         print("Interval:", U0, U1)
#         print(f)
#         print(df)
#         nonmonotone_funcs.append((f,(U0,U1)))
#         print()
#     return output

