# import fmodpy
# box_eval = fmodpy.fimport("box_eval.f90", verbose=True,
#                           module_link_args=["-lblas", "-llapack", "-lgfortran"])

# import numpy as np

# X = np.array([[1.,0.,1.],
#               [0.,1.,1.]], order="F")
# nu = np.array([2.,2.,2.], order="F")
# (xx,yy) = np.meshgrid( (np.arange(1,21)-2)/5,
#                        (np.arange(1,21)-2)/5 )
# p = np.asarray(np.vstack((yy.flatten(), xx.flatten())).T, order="F")

# print()
# print(X.shape)
# print(X)
# print()
# print(nu.shape)
# print(nu)
# print()
# print(p.shape)
# print()

# # # p = np.array(p[:5,:], order='F')
# p = np.array([[2.,2.]], order='F')
# b, error = box_eval.box_eval(X, nu, p)
# print("b: ")
# print("  ", b)
# exit()

# ====================================================================
# ====================================================================

import fmodpy
box_eval = fmodpy.fimport("box_eval.f90", verbose=True,
                          module_link_args=["-lblas", "-llapack", "-lgfortran"])

import numpy as np

# X = np.array([[1.,0.,1.],
#               [0.,1.,1.]], order="F")
# nu = np.array([1.,1.,1.], order="F") * 1

X = np.array([[1.,0., 1., 1.],
              [0.,1.,-1., 1.]], order="F")
nu = np.array([1.,1.,1.,1.], order="F") * 1


def f(x):
    p = np.asarray(x, order='F')
    b, e = box_eval.box_eval(X, nu, p)
    return b

from util.plotly import Plot
p = Plot()                             
padding = .02
# Get the x min and max
x_min_max = [sum(np.where(nu*X[0,:] < 0, nu*X[0,:], 0)), sum(np.where(nu*X[0,:] > 0, nu*X[0,:], 0))]
x_min_max[0] -= (x_min_max[1] - x_min_max[0]) * padding
x_min_max[1] += (x_min_max[1] - x_min_max[0]) * padding
# Get the y min and max
y_min_max = [sum(np.where(nu*X[1,:] < 0, nu*X[1,:], 0)), sum(np.where(nu*X[1,:] > 0, nu*X[1,:], 0))]
y_min_max[0] -= (y_min_max[1] - y_min_max[0]) * padding
y_min_max[1] += (y_min_max[1] - y_min_max[0]) * padding
# Create the plot
p.add_func("Box Spline", f, x_min_max, y_min_max, vectorized=True)
p.show(file_name="result.html")


# ====================================================================
# ====================================================================


# # print(b)
# # print(error)

# def f(x):
#     x = np.asarray(x[None,:], order='F')
#     b, _ = box_eval.box_eval(X, nu, x)
#     return b[0]

# from util.plotly import Plot
# p = Plot()                             
# p.add_func("Box Spline", f, *([[-.1,4.1]]*2))
# p.show(z_range=[-.2,.7], file_name="result.html")


#   #  5   2
#   # -0.20000  -0.20000
#   # -0.20000   0.00000
#   # -0.20000   0.20000
#   # -0.20000   0.40000
#   # -0.20000   0.60000
#   #  2   3
#   #  0.66667  -0.33333   0.33333
#   # -0.33333   0.66667   0.33333
