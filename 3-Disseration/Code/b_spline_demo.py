import numpy as np
import fmodpy
from util.plot import Plot

splines = fmodpy.fimport("splines.f08", verbose=True,
                         autocompile_extra_files=True,
                         module_link_args=["-lblas", "-llapack"])


unique_knots = [0., .7, 1., 2.]
min_knots = 2
max_knots = 7
eval_points = 5000

for num_knots in range(min_knots, max_knots+1):
    p = Plot(f"{'' if num_knots-3 < 0 else 'C'+str(num_knots-3)} B-Splines of order {num_knots-1}  (degree {num_knots-2}) over {unique_knots}")
    # Add the knots.
    for i,k in enumerate(unique_knots):
        p.add(f"Knot {i+1} at {k}", [k,k], [0,1], mode="lines", group="knots",
              color=p.color((0,0,0,.3)), dash="dot")
    print()
    print("order:  ",num_knots-1)
    print("degree: ",num_knots-2)
    ys = []
    x = np.linspace(min(unique_knots),max(unique_knots),eval_points)
    for start in range(2-num_knots, len(unique_knots)-1 ):
        # Build the list of knonts by repeating the first and last elements.
        knots = [unique_knots[max(0,min(len(unique_knots)-1,i))]
                 for i in range(start, start+num_knots)]
        print("start: ",start, knots)
        y = splines.eval_bspline(np.array(knots), x.copy())
        ys.append(y)
        # print("  ",x)
        # print("  ",y)
        p.add(str(knots), x, y, mode="lines")
    p.add("Sum", x, np.array(ys).sum(axis=0), mode="lines", 
          dash="dash", color=p.color((0,0,0)))
    # p.show(file_name="b-spline-demo.html")
    p.show(file_name="b-spline-demo.html",
           append=(num_knots>min_knots), 
           show=(num_knots==max_knots-1))

