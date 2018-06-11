import fmodpy
cbspline = fmodpy.fimport("cbspline.f90", verbose=True)
cbsplev = cbspline.cbsplev
cboxsplev = cbspline.cboxsplev

from util.plotly import Plot
import numpy as np

degree = 2

def bspline(x, degree, index=0, knots=None):
    # Initialize knots for a cardinal b-spline if none are given
    if (type(knots) == type(None)): knots = list(range((degree+1)+1))
    # Body of function
    if (degree == 0):
        # Base case, 0th degree, return 1 if in knot range else 0
        return int( knots[index] < x <= knots[index+1] )
    else:
        # Recursive case
        ramp_up = (x - knots[index]) / (knots[index+degree] - knots[index])
        ramp_down = (knots[1+index+degree] - x) / (knots[1+index+degree] - knots[1+index])
        return ( ramp_up * bspline(x, degree-1, index, knots) +
                 ramp_down * bspline(x, degree-1, index+1, knots) )


p = Plot()
# p.add_func("BSpline", lambda x: bspline(x, degree, knots=[0,.5,1.75,3]), [0,degree+1])
# p.add_func("BSpline", lambda x: bspline(x, degree), [0,degree+1])
p.add_func("Fort BSPline", lambda x: cbsplev(x, degree), [0,degree+1])
p.plot(file_name="spline.html")

# lower = np.zeros(2, order='F')
# upper = np.ones(2, order='F')
# def eval_box(x):
#     x = np.array(x, order='F')
#     return cboxsplev(x, degree, lower, upper)

# p = Plot()
# p.add_func("Box Spline", eval_box, [-.1,1.1], [-.1,1.1], plot_points=4000)
# p.plot(file_name="spline.html", append=True)
