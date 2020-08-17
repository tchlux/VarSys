# Custom exception for when there are data errors.
class BadUsage(Exception): pass
class BadData(Exception): pass
class BadValues(Exception): pass
class UnexpectedError(Exception): pass


# Given three points, this will solve the equation for the linear
# function which interpolates the first 2 values. Returns coefficient 2-tuple.
def linear(x, y):
    if len(x) != len(y): raise(BadUsage("X and Y must be the same length."))
    if len(x) < 2:      raise(BadUsage(f"At least 2 (x,y) coordinates must be given, received '{x}'."))
    x1, x2 = x[0], x[-1]
    y1, y2 = y[0], y[-1]
    slope = (y2 - y1) / (x2 - x1)
    return (y1 - slope*x1, slope)


# Given three points, this will solve the equation for the quadratic
# function which interpolates first 3 values. Returns coefficient 3-tuple.
def quadratic(x, y):
    if len(x) != len(y): raise(BadUsage("X and Y must be the same length."))
    if len(x) < 3:      raise(BadUsage(f"At least 3 (x,y) coordinates must be given, received '{x}'."))
    x1, x2, x3 = x[0], x[len(x)//2], x[-1]
    y1, y2, y3 = y[0], y[len(y)//2], y[-1]
    a = -((-x2 * y1 + x3 * y1 + x1 * y2 - x3 * y2 - x1 * y3 + x2 * y3)/((-x1 + x2) * (x2 - x3) * (-x1 + x3)))
    b = -(( x2**2 * y1 - x3**2 * y1 - x1**2 * y2 + x3**2 * y2 + x1**2 * y3 - x2**2 * y3)/((x1 - x2) * (x1 - x3) * (x2 - x3)))
    c = -((-x2**2 * x3 * y1 + x2 * x3**2 * y1 + x1**2 * x3 * y2 - x1 * x3**2 * y2 - x1**2 * x2 * y3 + x1 * x2**2 * y3)/((x1 - x2) * (x1 - x3) * (x2 - x3)))
    return (c,b,a)


# Given data, determine a polynomial fit to that data.
def derivatives(x, y, left=None, right=None, method=linear):
    # Set the left and right derivatives if they are not provided.
    if (type(left) == type(None)):  left  = [0]*(len(x)-1)
    if (type(right) == type(None)): right = [0]*(len(x)-1)
    # Extend the left and right derivatives if they are not long enough.
    if (len(left) < len(x)):  left  = left  + [0]*(len(x)-1 - len(left))
    if (len(right) < len(x)): right = right + [0]*(len(x)-1 - len(right))
    # Compute all other derivatives with quadratic fits.
    derivatives = [left] + [ list() for i in range(len(x) - 2) ] + [right]
    for d in range(len(x)-1):
        print("d: ",d)
        for i in range(1, len(x)-1):
            print(" i: ", i)
            positions = (x[i-1], x[i], x[i+1])
            if (d == 0): values = (y[i-1], y[i], y[i+1])
            else:        values = (derivatives[i-1][d-1],
                                   derivatives[i][d-1],
                                   derivatives[i+1][d-1])
            # Compute the derivative, get the linear term of the fit.
            derivatives[i].append( method(positions, values)[1] )
    return derivatives


x = [0,1,2]
y = [1,2,3]
out = derivatives(x, y, method=quadratic)

import numpy as np

print()
print(np.array(out))

from util.plot import Plot
f = lambda x: x**3 - 3*x**2 + 3*x + 1
p = Plot()
p.add("Points", x, y)
p.add_func("Fit", f, [-1, 3])
p.show()
