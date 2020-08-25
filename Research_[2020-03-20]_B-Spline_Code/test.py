import numpy as np
import og_fmodpy as fmodpy
from util.plot import Plot

splines = fmodpy.fimport("splines.f90", verbose=True,
                         autocompile_extra_files=True,
                         module_link_args=["-lblas", "-llapack"])


TEST_BANDED_MATRIX = False
TEST_B_SPLINE = True
TEST_SPLINE = False


# =============================================
#      Test the packing of a banded matrix     
# =============================================

if TEST_BANDED_MATRIX:
    # AB(KL+KU+1+i-j,j) = A(i,j) for max(1,j-KU)<=i<=min(N,j+KL)

    N = 6
    KL = 3
    KU = 3
    NRHS = 1

    LDAB = 2*KL + KU + 1
    AB = np.zeros((LDAB,N))
    IPIV = np.zeros(N)
    LDB = N
    B = -(np.arange(LDB)+1)

    A = np.zeros((N,N))
    for i in range(1,A.shape[0]+1):
        for j in range(1,A.shape[1]+1):
            if (max(1,j-KU) <= i <= min(N,j+KL)):
                A[i-1,j-1] = i*10 + j
                AB[KL+KU+1+i-j-1,j-1] = i*10 + j

    print()
    print(A)
    print()
    print(AB)
    print()


# ===========================================
#      Test the B-spline evaluation code     
# ===========================================

if TEST_B_SPLINE:
    # knots = np.array([0.,1,3,3,6])+1
    # knots = np.array([1.,1,1,1])
    # knots = np.array([0.,0,0,1])
    # knots = np.array([0.,1,1,1]) 
    # knots = np.array([0.,0,0,1,1])
    # knots = np.array([0.,1,1,2]) + 1
    # knots = np.array([0.,0,0,1,1,1,2,2,2])
    r = 2
    knots = np.array([0.] + [1]*r + r*[2])
    knots *= len(knots)-2

    knots = np.array([0.0000000000000000,        0.0000000000000000,        0.0000000000000000,        0.0000000000000000,        0.0000000000000000,        0.0000000000000000,       0.14285714285714285])


    print()
    print("len(knots) : ",len(knots))
    print("      knots: ",knots)

    # Generate evaluation points that hit every single knot.
    x = [knots[0]-.1]
    for i in range(len(knots)-1):
        if (knots[i] == knots[i+1]): continue
        x += list(np.linspace(knots[i], knots[i+1], 20))[:-1]
    x = np.array(x + [knots[-1], knots[-1]+.1], dtype=float)
    # Evaluate the spline value everywhere.
    y, info = splines.eval_bspline(knots, x.copy(), d=0)
    print("info: ",info)

    # Use exact-arithmetic to compute the "truth" derivatives and integrals.
    from polynomial import BSpline
    from fraction import Fraction
    truth = BSpline( list(map(Fraction,knots)) )

    # Compute values at knots.
    round_to = 3
    max_d_vals = []
    for d in range(len(knots)):
        at_knots, info = splines.eval_bspline(knots,knots.copy(), d=d)
        at_true  = np.array([truth(v, d=d) for v in knots])
        if d != 0: max_d_vals.append(max(abs(at_knots)))
        print()
        print(f"D{d} at knots: ", " ".join([f"{v: .{round_to}f}" for v in at_knots]))
        print(f"      truth: ", " ".join([f"{v: .{round_to}f}" for v in at_true]))

    for d in range(1,1): #len(knots)):
        at_knots, info = splines.eval_bspline(knots,knots.copy(), d=-d)
        at_true  = np.array([truth(v, d=-d) for v in knots])
        print()
        print(f"I{d} at knots: ", " ".join([f"{v: .{round_to}f}" for v in at_knots]))
        print(f"      truth: ", " ".join([f"{v: .{round_to}f}" for v in at_true]))
    print()

    # Compute the maximum derivative values observed.
    max_d_vals = np.array(max_d_vals, dtype=float)
    to_show = np.argsort(max_d_vals)
    # print("max_d_vals          : ",max_d_vals)
    # print("to_show             : ",to_show)
    # print("max_d_vals[to_show] : ",max_d_vals[to_show])


    # Make a pretty picture of the B-spline and its derivatives.
    styles = ["dot", "dash", "dashdot"]
    # Add the B-spline function itself.
    continuity = "" if (len(knots) <= 2) else f"C{len(knots)-3} "
    p = Plot(continuity+f"B-Spline, piecewise order {len(knots)-1} polynomial")
    p.add("B-spline", x, y, mode="markers+lines", color=p.color(1), group='truth')
    p.add("Truth", x, [truth(v) for v in x], mode="markers+lines",
          color=p.color(1, alpha=.2), fill="toprevy", group='truth')

    # Add all interesting derivatives.
    # for d in range(max(1,len(knots)-3), len(knots)):
    for d in to_show[-4:]:
        d = d + 1
        color = p.color(d-1) if d == 1 else p.color(d)
        dash = styles[(d-1)%len(styles)]
        dy, info = splines.eval_bspline(knots,x.copy(), d=d)
        name = f"D{d}"
        p.add(name, x, dy, mode="lines", color=color, dash=dash, group=name)
        # Add the true derivataive.
        color = p.color(d-1,alpha=.3) if d == 1 else p.color(d,alpha=.3)
        p.add(f"True {name}", x, [truth(v,d=d) for v in x], color=color,
              group=name, mode="lines", fill="toprevy")

    # Add all interesting integrals.
    for d in range(max(1,len(knots)-3), 1): #len(knots)):
        i = d + len(knots)-2
        color = p.color(i)
        dash = styles[(i-1)%len(styles)]
        dy, info = splines.eval_bspline(knots,x.copy(), d=-d)
        name = f"I{d}"
        p.add(name, x, dy, mode="lines", color=color, dash=dash, group=name)
        # Add the true integral.
        color = p.color(i, alpha=.3)
        p.add(f"True {name}", x, [truth(v,d=-d) for v in x], color=color,
              group=name, mode="lines", fill="toprevy")

    # Add the knots themselves.
    p.add("Knot Locations", [None], [None], color='rgba(0,0,0,.1)', group='k')
    for i in range(len(knots)):
        p.add("", [knots[i],knots[i]], [-1,1],
              color='rgba(0,0,0,.1)', group='k',
              show_in_legend=False, mode="lines")

    # Generate the visual.
    p.show(file_name=f"deriv_test-{len(knots)-1}.html") #, y_range=[-3*max(y),4*max(y)])


# ==================================
#      Test the spline-fit code     
# ==================================

# import splines
# from util.plot import Plot
# knots = [0.5488135 , 0.5488135 , 0.5488135 , 0.5488135 , 0.5488135 ,
#        0.5488135 , 0.5488135 , 0.5488135 , 0.5488135 , 0.5488135 ,
#        0.5488135 , 0.5488135 , 0.5488135 , 0.5488135 , 0.5488135 ,
#        0.5488135 , 0.5488135 , 0.5488135 , 0.5488135 , 0.5488135 ,
#        0.71518937, 0.71518937, 0.71518937, 0.71518937, 0.71518937,
#        0.71518937, 0.71518937, 0.71518937, 0.71518937, 0.71518937,
#        0.71518937, 0.71518937, 0.71518937, 0.71518937, 0.71518937,
#        0.71518937, 0.71518937, 0.71518937, 0.71518937, 0.71518937]

# coefs = [0.76604712, 0.76975691, 0.77350211, 0.7772835 , 0.78110183,
#        0.7849579 , 0.78885247, 0.79278634, 0.79676031, 0.80077518,
#        0.49655397, 0.50164066, 0.50679749, 0.51202474, 0.51732269,
#        0.52269164, 0.52813186, 0.53364366, 0.53922733, 0.54488318]

# for d in range( (len(knots)-len(coefs))+1 ):
#     p = Plot(f"d{d}")
#     for i in range(len(coefs)):
#         k = np.array(knots[i:i+(len(knots)-len(coefs))+1], dtype=float)
#         x = np.linspace(knots[0], knots[-1], 1000)
#         y = splines.eval_bspline(k, x.copy(), d=d)
#         p.add(f"B{i+1}", x, y, mode="lines")
#     p.show(append=True)
# exit()



if TEST_SPLINE:
    RANDOMIZE_TEST = True
    VISUALIZE_TEST = True
    # 
    # TODO: Large multipliers on the knot spacing cause the numerical
    #       stability to favor the higher derivatives.
    # 
    #       Medium multipliers on the knot spacing cause the stability
    #       to favor the value (and lower derivatives).
    # 
    #       Small multipliers on the knot spacing cause the stability
    #       to only be good for the function values.
    # 
    #       This must have to do with the relative size of the
    #       derivatives to the function values, is there a way for me
    #       to achieve the "balanced" stability in all outcomes? 
    #       I will have to find a way to make the magnitudes of the
    #       derivative evaluations and function evaluations similar
    #       for the underlying B-splines, there must be a way to do
    #       that.
    # 
    #       I can find multipliers (no internal rescaling) that cause
    #       the accuracy to be "level" across all derivatives. It
    #       does not seem to involve the ratio between the B-spline
    #       coefficients. I suspect it involves the coefficient of
    #       difference between the function and derivative values.
    # 
    #       Linearly rescaaling the values appears to have no
    #       difference between the (order of magnitude) gap between
    #       the largest -> smallest error in the approximation of all
    #       derivatives. Yet, rescaling the knots does appear to have
    #       an effect. The nonlinear relationship between knot scaling
    #       and derivative values must be important.
    # 
    #       Rescaling the knot spacing to be normalized in the code
    #       has no effect on the stability. Rescaling the knots before
    #       calling the code has an effect on stability.
    # 
    #       The only thing that consistently works is scaling up knots
    #       without changing values (changing the approximation
    #       problem). Once the MAX(A) is <= 1.0, the accuracy starts
    #       to flip. The higher derivatives become more accurate than
    #       the function value. Eventually the function value in
    #       unacceptably inaccurate.
    # 
    if RANDOMIZE_TEST:
        np.random.seed(0)
        num_knots = 3
        continuity = 1
        # Generate some random knots and values.
        knots = np.random.random(size=(num_knots))
        knots.sort()
        values = 10 * (2*(np.random.random(size=(num_knots, continuity+1)) - 1/2))
        # Normalize the knots and values differently.
        # knots *= 100
        # knots -= min(knots)
        # knots /= max(knots)
        # knots /= min(knots[1:] - knots[:-1])
        # knots *= 10
        # knots /= min(knots[1:] - knots[:-1])
        # values /= 100000
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

