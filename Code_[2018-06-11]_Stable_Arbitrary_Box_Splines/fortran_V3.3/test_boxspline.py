import numpy as np
import fmodpy
boxspline = fmodpy.fimport("boxspline.f90", verbose=True,
    module_link_args=["-lblas", "-llapack", "-lgfortran"])


for PTS in ["1", "2K", "4K"]:
    for ELEM in ["TenP", "Diag", "ZP"]:
        for mult in [1,2,3]:
            name = f"_{PTS}_{ELEM}_{mult}.csv"
            print(f"Starting {name}", end="  ...   ")
            # ====================================================================
            if ELEM == "TenP":
                # Bi-linear / bi-quadratic / bi-cubic function
                dvecs = np.array([[1.,0.],
                                  [0.,1.]], order="F")
                mults = np.array([1 ,1 ], order="F", dtype=np.int32) * mult
            if ELEM == "Diag":
                # Not sure if this is named, but it's another test.
                dvecs = np.array([[1.,0.,1.],
                                  [0.,1.,1.]], order="F")
                mults = np.array([1 ,1 ,1 ], order="F", dtype=np.int32) * mult
            if ELEM == "ZP":
                # Zwart-Powell element
                dvecs = np.array([[1.,0., 1., 1.],
                                  [0.,1.,-1., 1.]], order="F")
                mults = np.array([1 ,1 ,1 ,1 ], order="F", dtype=np.int32) * mult
            # ====================================================================

            if PTS == "1":
                eval_pts = np.asarray(np.asarray([np.sum(dvecs*mult,axis=1)/2]).T, order='F', dtype=np.float64)
            if PTS == "2K":
                eval_pts = np.meshgrid(list(range(50)), list(range(50)))
                eval_pts = np.vstack((eval_pts[0].flatten(), eval_pts[1].flatten())).T
                eval_pts = np.asarray(eval_pts.T, order='F', dtype=np.float64)
            if PTS == "4K":
                def f(x):
                    global eval_pts
                    eval_pts = np.asarray(x.T, order='F')
                    return np.sum(eval_pts, axis=0)

                from util.plot import Plot
                p = Plot()                             
                padding = .05
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
                p.add_func("Box Spline", f, x_min_max, y_min_max, vectorized=True,
                           use_gradient=True, plot_points=4000)

            # ====================================================================

            np.savetxt(f"dvecs"    +name, dvecs,    delimiter=",", header=",".join(map(str,dvecs.shape)),    comments="")
            np.savetxt(f"mults"    +name, mults,    delimiter=",", header=",".join(map(str,mults.shape)),    comments="")
            np.savetxt(f"eval_pts" +name, eval_pts, delimiter=",", header=",".join(map(str,eval_pts.shape)), comments="")
            print("done.")
