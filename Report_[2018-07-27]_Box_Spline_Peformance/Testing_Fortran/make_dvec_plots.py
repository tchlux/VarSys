import numpy as np
from util.plot import Plot, multiplot
import boxspline

tp = np.array([[1.,0.],
               [0.,1.]], order="F")
cr = np.array([[1.,0.,1.],
               [0.,1.,1.]], order="F")
zp = np.array([[1.,0., 1.,1.],
               [0.,1.,-1.,1.]], order="F")

nms, els, mls = ("TenP","Cour","ZP"), (tp, cr, zp), ((1,2,3), (1,2), (1,2))
for name, dvec, mults in zip(nms, els, mls):
    plots = []
    for mult in mults:
        p = Plot(f"","","","")
        mins = np.sum(np.where(dvec < 0, dvec, 0), axis=1) * mult
        maxs = np.sum(np.where(dvec > 0, dvec, 0), axis=1) * mult
        mins -= .15 * (maxs - mins)
        maxs += .15 * (maxs - mins)
        mult = np.ones(dvec.shape[1], dtype=np.int32) * mult
        # Define a function that evaluates the box-spline.
        def box(x):
            x = np.array(x, order="F")
            evals, error = boxspline.boxsplev(dvec, mult, x)
            if (error != 0): print("ERROR:", error)
            return evals
        p.add_func("", box, *list(zip(mins,maxs)), plot_points=4000,
                   show_in_legend=False, vectorized=True, use_gradient=True)
        plots.append( p.plot(z_range=[-.05,2.05], show=False, html=False) )
    multiplot([plots], file_name=f"{name}.html", show=False)
