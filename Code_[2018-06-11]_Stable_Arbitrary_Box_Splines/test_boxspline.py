import numpy as np
import fmodpy
boxspline = fmodpy.fimport("boxspline.f90", verbose=True,
    module_link_args=["-lblas", "-llapack", "-lgfortran"])

# dvecs = np.array([[1.,0.,1.],
#                   [0.,1.,1.]], order="F")
# mults = np.array([1 ,1 ,1 ], order="F", dtype=np.int32) * 2

dvecs = np.array([[1.,0., 1., 1.],
                  [0.,1.,-1., 1.]], order="F")
mults = np.array([1 ,1 ,1 ,1 ], order="F", dtype=np.int32) * 1

def f(x):
    eval_pts = np.asarray(x, order='F')
    box_evals, error, info = boxspline.boxsplev(dvecs, mults, eval_pts)
    if (error != 0):
        print()
        print( "Error evaluating box spline:")
        print(f"   ERROR = {error:02d}")
        print(f"   INFO  = {info:02d}")        
        print()
    return box_evals

from util.plotly import Plot
p = Plot()                             
padding = .02
# Get the x min and max
x_min_max = [sum(np.where(mults*dvecs[0,:] < 0, mults*dvecs[0,:], 0)), 
             sum(np.where(mults*dvecs[0,:] > 0, mults*dvecs[0,:], 0))]
x_min_max[0] -= (x_min_max[1] - x_min_max[0]) * padding
x_min_max[1] += (x_min_max[1] - x_min_max[0]) * padding
# Get the y min and max
y_min_max = [sum(np.where(mults*dvecs[1,:] < 0, mults*dvecs[1,:], 0)), 
             sum(np.where(mults*dvecs[1,:] > 0, mults*dvecs[1,:], 0))]
y_min_max[0] -= (y_min_max[1] - y_min_max[0]) * padding
y_min_max[1] += (y_min_max[1] - y_min_max[0]) * padding
# Create the plot
p.add_func("Box Spline", f, x_min_max, y_min_max, vectorized=True)
p.show(file_name="result.html")

# C-x C-k e -- Edit currently defined keyboard macro.
# C-x C-k n -- Give name to the most recently defined keyboard macro (session).
# C-x C-k b -- Bind most recent keyboard macro to a key (session).

