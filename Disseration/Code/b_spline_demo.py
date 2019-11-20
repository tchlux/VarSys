import numpy as np
import fmodpy
from util.plot import Plot

splines = fmodpy.fimport("splines.f08", verbose=True,
                         autocompile_extra_files=True,
                         module_link_args=["-lblas", "-llapack"])


for size in range(2,5+1):
    p = Plot(f"B-Splines of order {size-1}  (degree {size-2})")
    print()
    print("size: ",size)
    for start in range(1,size):
        knots = np.array([0.] * (size-start) + [1.] * start)
        x = np.linspace(0,1,100)
        y = splines.eval_bspline(knots, x.copy())
        print(" knots: ",knots)
        # print("  ",x)
        # print("  ",y)
        p.add(str(start), x, y, mode="lines")
    # p.show(file_name="b-spline-demo.html")
    p.show(file_name="b-spline-demo.html", append=(size>2), show=(size==5))

