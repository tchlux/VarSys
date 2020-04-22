
import fmodpy
splines = fmodpy.fimport("splines.f90", verbose=True,
                         autocompile_extra_files=True,
                         module_link_args=["-lblas", "-llapack"])


TEST_B_SPLINE = False
TEST_SPLINE = True


# ==================================
#      Test the B-spline code     
# 
if TEST_B_SPLINE:
    import numpy as np
    from util.plot import Plot
    left, right = 0.5488135, 0.71518937
    knots = [left]*20 + [right]*20
    coefs = [0.76604712, 0.76975691, 0.77350211, 0.7772835 , 0.78110183,
             0.7849579 , 0.78885247, 0.79278634, 0.79676031, 0.80077518,
             0.49655397, 0.50164066, 0.50679749, 0.51202474, 0.51732269,
             0.52269164, 0.52813186, 0.53364366, 0.53922733, 0.54488318]

    for d in range( (len(knots)-len(coefs))+1 ):
        p = Plot(f"d{d}")
        for i in range(len(coefs)):
            k = np.array(knots[i:i+(len(knots)-len(coefs))+1], dtype=float)
            x = np.linspace(knots[0], knots[-1], 1000)
            y, status = splines.eval_bspline(k, x.copy(), d=d)
            p.add(f"B{i+1}", x, y, mode="lines")
        first = (d == 0)
        last = (d == (len(knots)-len(coefs)))
        p.show(append=True, show=first or last)


# ==================================
#      Test the spline-fit code     
# 
if TEST_SPLINE:
    import numpy as np
    RANDOMIZE_TEST = True
    VISUALIZE_TEST = True
    if RANDOMIZE_TEST:
        np.random.seed(1)
        num_knots = 5
        continuity = 2
        # Generate some random knots and values.
        knots = 10 * np.random.random(size=(num_knots))
        knots.sort()
        values = 10 * (2*(np.random.random(size=(num_knots, continuity+1)) - 1/2))
        print()
        print("knots:")
        print(knots)
        print()
        print("values:")
        print(values)
        print()
        # Get the spline knots and spline coefficients.
        values = np.asfortranarray(values)
        sk, sc, info = splines.fit_spline(knots, values)
        print("info: ",info)
        print()
        # print("order: ",len(sk) - len(sc) + 1)
        # print()
        # print("sk:")
        # print(repr(sk))
        print("sc:")
        print(repr(sc))
        print()
        # sc *= -1
        for d in range(continuity+1):
            y, info = splines.eval_spline(sk, sc, knots.copy(), d=d)
            error = abs(y - values[:,d])
            print(f"f^({d:2d}) max error:   {max(error):.3e}   at  ", knots[error>=max(error)*.9])
        print()

        if VISUALIZE_TEST:
            from util.plot import Plot
            # Make a pretty picture of the B-spline and its derivatives.
            padding = 0 #.1
            x = np.linspace(min(knots)-padding,max(knots)+padding,max(1000,(len(knots)-1)*10+1))
            y, info = splines.eval_spline(sk, sc, x.copy(), d=0)
            p = Plot("Polynomial Interpolant")
            p.add("Spline", x, y, mode="lines", color=p.color(1))
            # --------------------------------------------------------
            from polynomial import Spline
            from fraction import Fraction
            def exact(v):
                try:    return list(map(exact,v))
                except: return Fraction(v)
            truth = Spline(exact(knots), exact(values))
            true_y = list(map(truth, x))
            p.add("Truth", x, true_y, mode="lines", color=p.color(1), dash="dot",fill="toprevy")
            # --------------------------------------------------------
            p.add("Knots", knots, values[:,0], color=p.color(1))
            p.show(file_name=f"spline_test-N{len(values)}-C{len(values[0])}.html")

    else:
        # values = np.array([[1.,0,2,1],
        #                    [.5,0,1,2],
        #                    [0,0,0,3]]).T
        values = np.array([[1., 0],
                           [0., 0],
                           [1., 0],
                           [0., 0]])
        knots = np.linspace(0,values.shape[0]-1, values.shape[0])
        continuity = values.shape[1]
        print("knots:  ")
        print(knots)
        print("values:")
        print(values)
        # Get the spline knots and spline coefficients.
        values = np.asfortranarray(values)
        sk, sc, info = splines.fit_spline(knots, values)

        # sc *= -1
        for d in range(continuity):
            y, info = splines.eval_spline(sk, sc, knots.copy(), d=d)
            error = abs(y - values[:,d])
            print(f"f^({d}) max error: {max(error):.3e}")
        print()
        print("sk: ",sk)
        print("sc: ",sc)
        print()
        print("Knots:")
        print("   ",knots)
        for d in range(values.shape[1]*2):
            print(f"{d}: ",splines.eval_spline(sk, sc, knots.copy(), d=d)[0])
        print()

        if VISUALIZE_TEST:
            from util.plot import Plot
            padding = 0 #.1
            x = np.linspace(min(knots)-padding,max(knots)+padding,1000)
            print()
            print(f"evaluating spline..")
            y, info = splines.eval_spline(sk, sc, x.copy(), d=0)
            # Make a pretty picture of the B-spline and its derivatives.
            p = Plot("Hermite Interpolant")
            p.add("Spline", x, y, mode="lines", color=p.color(1), group="s")
            # --------------------------------------------------------
            from polynomial import Spline
            from fraction import Fraction
            def exact(v):
                try:    return list(map(exact,v))
                except: return Fraction(v)
            truth = Spline(exact(knots), exact(values))
            true_y = list(map(truth, x))
            p.add("Truth", x, true_y, mode="lines", color=p.color(1), dash="dot",fill="toprevy")
            # --------------------------------------------------------
            p.add("Knots", knots, values[:,0], color=p.color(1), group="s")
            # Add all interesting derivatives.
            styles = ["dot", "dash", "dashdot"]
            for d in range(1, values.shape[1]*2):
                color = p.color(d-1) if d == 1 else p.color(d)
                dash = styles[(d-1)%len(styles)]
                print()
                print(f"{d} derivative evaluating..")
                dy, info = splines.eval_spline(sk, sc, x.copy(), d=d)
                p.add(f"D{d}", x, dy, mode="lines", color=color, dash=dash,group=d)
                if (d < values.shape[1]):
                    p.add(f"k{d}", knots, values[d,:], color=color, group=d)
            y, info = splines.eval_spline(sk, sc, x.copy(), d=-1)
            p.add("Integral", x, y, mode="lines", color=p.color(0), group="i")

            p.show(file_name=f"spline_test-N{len(values)}-C{len(values[0])}.html")

