import numpy as np
import fmodpy
boxspline = fmodpy.fimport("boxspline.f90", verbose=True,
    module_link_args=["-lblas", "-llapack", "-lgfortran"])

# dvecs = np.array([[1.,0.],
#                   [0.,1.]], order="F")
# mults = np.array([1 ,1 ], order="F", dtype=np.int32) * 2

# dvecs = np.array([[1.,0.,1.],
#                   [0.,1.,1.]], order="F")
# mults = np.array([1 ,1 ,1 ], order="F", dtype=np.int32) * 2

dvecs = np.array([[1.,0., 1., 1.],
                  [0.,1.,-1., 1.]], order="F")
mults = np.array([1 ,1 ,1 ,1 ], order="F", dtype=np.int32) * 1

# ====================================================================

eval_pts = np.meshgrid(list(range(50)), list(range(50)))
eval_pts = np.vstack((eval_pts[0].flatten(), eval_pts[1].flatten())).T
eval_pts = np.asarray(eval_pts[:1], order='F', dtype=np.float64)

import time
start = time.time()
box_evals, error = boxspline.boxsplev(dvecs, mults, eval_pts)
total = time.time() - start
print(box_evals, error)
print(f" {total:.2e} second evaluation time at {eval_pts.shape[0]} points..")
print()
print((58.4973 - 0.00644803) / total) # Matlab execution time minus allocation time.
exit()

# ====================================================================

# def f(x):
#     eval_pts = np.asarray(x, order='F')
#     box_evals, error = boxspline.boxsplev(dvecs, mults, eval_pts)
#     if (error != 0):
#         print()
#         print( "Error evaluating box spline:")
#         print(f"   ERROR = {error:02d}")
#         print()
#     return box_evals

# from util.plot import Plot
# p = Plot()                             
# padding = .05
# # Get the x min and max
# x_min_max = [sum(np.where(mults*dvecs[0,:] < 0, mults*dvecs[0,:], 0)), 
#              sum(np.where(mults*dvecs[0,:] > 0, mults*dvecs[0,:], 0))]
# x_min_max[0] -= (x_min_max[1] - x_min_max[0]) * padding
# x_min_max[1] += (x_min_max[1] - x_min_max[0]) * padding
# # Get the y min and max
# y_min_max = [sum(np.where(mults*dvecs[1,:] < 0, mults*dvecs[1,:], 0)), 
#              sum(np.where(mults*dvecs[1,:] > 0, mults*dvecs[1,:], 0))]
# y_min_max[0] -= (y_min_max[1] - y_min_max[0]) * padding
# y_min_max[1] += (y_min_max[1] - y_min_max[0]) * padding
# # Create the plot
# p.add_func("Box Spline", f, x_min_max, y_min_max, vectorized=True,
#            use_gradient=True, plot_points=4000)
# p.show(z_range=[-.05,2.05], file_name="result.html")

# C-x C-k e -- Edit currently defined keyboard macro.
# C-x C-k n -- Give name to the most recently defined keyboard macro (session).
# C-x C-k b -- Bind most recent keyboard macro to a key (session).

