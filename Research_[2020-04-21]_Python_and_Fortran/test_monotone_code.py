import numpy as np
from polynomial import Spline, fit, local_polynomial, local_quadratic, \
    polynomial, Polynomial
# from monotone import monotone_quintic_spline
from search_monotone import monotone_quintic_spline
from fraction import Fraction
from util.plot import Plot



# Make a plot of the points where the constraints 
def example_points(x, y, ms=10, p=None):
    # Find all "zero derivative" points.
    d_pos = []
    d_neg = []
    d_zero = []
    dd_zero = []
    for i in range(len(x)):
        left  = (y[i]-y[i-1]) if (i > 0)        else ((-1)**int(y[0] > y[1])   if len(y) > 1 else 0)
        right = (y[i+1]-y[i]) if (i+1 < len(y)) else ((-1)**int(y[-2] > y[-1]) if len(y) > 1 else 0)
        if   (left * right) < 0:             d_zero.append(i)
        elif ((left == 0) or (right == 0)):  dd_zero.append(i)
        elif (left > 0):                     d_pos.append(i)
        elif (left < 0):                     d_neg.append(i)
        else: print("This is unexpected..", left, right)
    # Create the plot.
    if (p is None): p = Plot()
    p.add("Local Extrema (Df = 0)", [x[i] for i in d_zero], [y[i] for i in d_zero], color=0, symbol="circle", marker_size=ms, marker_line_width=2)
    p.add("Flat slope (DDf = 0, Df=0)", [x[i] for i in dd_zero], [y[i] for i in dd_zero], color=1, symbol="square", marker_size=ms, marker_line_width=2)
    p.add("Monotone increasing (Df ≥ 0)", [x[i] for i in d_pos], [y[i] for i in d_pos], color=2, symbol="triangle-up", marker_size=ms, marker_line_width=2)
    p.add("Monotone decreasing (Df ≤ 0)", [x[i] for i in d_neg], [y[i] for i in d_neg], color=3, symbol="triangle-down", marker_size=ms, marker_line_width=2)
    return p


# Visualize multiple methods.
def visualize_multiple(x, y, name=["Local quintic"], title="",
                       cubics=[False], color=[1], dash=[None],
                       **kwargs):
    # Initialize the plot.
    p = Plot(title=title, font_family="times")
    # Add visuals of all methods.
    for (n,is_cubic,c,d) in zip(name, cubics, color, dash):
        if not is_cubic:
            f = monotone_quintic_spline(x,y)
        else:
            from scipy.interpolate import PchipInterpolator
            f = PchipInterpolator(x, y)
        p.add_func(n, f, [min(x), max(x)], color=c, dash=d)
        # p.add_func("D "+n, f.derivative(), [min(x), max(x)], color=c, opacity=.2)

    # Add the points styled according to constraints (last to be on top).
    example_points(x,y, ms=8, p=p)
    # Generate the plot and return.
    min_max = [float(min(y)), float(max(y))]
    width = min_max[1] - min_max[0]
    return p.show(y_range=[min_max[0]-width/8, min_max[1]+width/8], **kwargs)



# Visualize the application of one approximation method.
def visualize_one(x, y, append=False, title=None, show=False,
                  file_name="pmqsi_example.html", local=True,
                  legend=False, cubic=False, color=1,
                  pure=False, verbose=False, p=None, dash=None):
    if (p is None): p = Plot(font_family="times")

    # Convert "x" and "y" to Fraction types for higher accuracy.
    x = list(map(Fraction, x))
    y = list(map(Fraction, y))

    # Add the cubic interpolant.
    if cubic:
        if dash is None:
            if cubic: file_name, dash = "pchip_example.html", None
            else:     dash = "dash"
        from scipy.interpolate import PchipInterpolator
        f = PchipInterpolator(x, y)
        p.add_func("Cubic Interpolant", f, [min(x), max(x)], dash=dash, color=2)
    # Add the quintic interpolant.
    else:
        # Create the spline.
        local_fits = {}
        f = monotone_quintic_spline(x,y, local_fits=local_fits,
                                    verbose=verbose)
        p.add_func("Quintic Interpolant", f, [min(x), max(x)],
                   color=color, dash=dash, plot_points=5000)
        # Add the local interpolants used to generate derivatives.
        if local:
            for i in range(len(x)):
                if i not in local_fits:
                    print(f"Local fits missing {i}, skipping..")
                    continue
                f, bounds = local_fits[i]
                if (i > 0): bounds[0] = (x[i-1] + x[i]) / 2
                if (i+1 < len(x)):  bounds[1] = (x[i+1] + x[i]) / 2
                p.add_func(f" local approx about {i}", f, bounds, dash="dot", opacity=0.4)

    # Add a pure "fit" that maintains desired continuity.
    if pure: p.add_func("pure quintic fit", fit(x, y, continuity=2), [min(x), max(x)], dash="dot")

    # Add the points styled according to constraints (last to be on top).
    example_points(x,y, ms=6, p=p)

    # "Piecwise monotone quintic spline interpolant<br>"+
    if (title is not None): p.title = title.title()

    # Generate the plot and return.
    min_max = [float(min(y)), float(max(y))]
    width = min_max[1] - min_max[0]
    return p.show(y_range=[min_max[0]-width/8, min_max[1]+width/8],
                  file_name=file_name, append=append,
                  show_legend=legend, show=show)




if True:

    cubics = [False] #, True]
    names = ["Quintic Facet", "PCHIP"]
    # colors = [(200,20,80), 3]
    colors = [1, 3]
    dashes = [None, None] #"dot"]
    legend = False
    show = True

    test = -1
    # --------------------------------------------------------------------
    # Watson test.
    y = [0,1,1,1,0,20,19,18,17,0,0,3,0,1,6,16,16.1,1]
    x = list(range(len(y)))

    test += 1
    visualize_multiple(x, y, name=names, cubics=cubics,
                       color=colors, dash=dashes,
                       show_legend=legend, show=show,
                       title=f"Test {test}",
                       file_name=f"test_{test}.html")


    # --------------------------------------------------------------------
    # Sin function.
    np.random.seed(1)
    x = np.random.random(size=(30,)) * 6*np.pi
    x.sort()
    y = np.sin(x)

    test += 1
    visualize_multiple(x, y, name=names, cubics=cubics,
                       color=colors, dash=dashes,
                       show_legend=legend, show=show,
                       title=f"Test {test}",
                       file_name=f"test_{test}.html")


    # --------------------------------------------------------------------
    # Large tangent test.
    large_tangent = lambda x: -(1 + ( 1/ (x-1) ))
    np.random.seed(1)
    x = np.linspace(0, .98, 10)
    y = [large_tangent(v) for v in x]
    test += 1
    visualize_multiple(x, y, name=names, cubics=cubics,
                       color=colors, dash=dashes,
                       show_legend=legend, show=show,
                       title=f"Test {test}",
                       file_name=f"test_{test}.html")
    

    # --------------------------------------------------------------------
    # Random monotone data test.
    np.random.seed(2)
    N = 16
    x = np.random.random(size=(N,))
    y = np.random.random(size=(N,))
    x.sort(); y.sort()
    test += 1
    visualize_multiple(x, y, name=names, cubics=cubics,
                       color=colors, dash=dashes,
                       show_legend=legend, show=show,
                       title=f"Test {test}",
                       file_name=f"test_{test}.html")


if False:
    # --------------------------------------------------------------------
    # Generate the first visuals explicitly requested by Dr. Watson.
    y = [0,1,1,1,0,20,19,18,0,0,1,0,1]
    y = [0,1,1,1,0,12,10,8,0,0,1,0,1]
    x = list(range(len(y)))

    # x = x[2:-3]
    # y = y[2:-3]

    # # --------------------------------------------------------------------
    # # Linear function with a small plateau.
    # small = 2**(-26)
    # x = [-small, 0, 1, 1.9, 2, 2.1, 3, 4, 4+small]
    # y = [-small, 0, 1, 1.99, 2, 2.01, 3, 4, 4+small]
    # x = x[1:-1]
    # y = y[1:-1]

    # # --------------------------------------------------------------------
    # # Absolute value function.
    # x = np.linspace(-1,1,11)
    # y = abs(x)

    # # --------------------------------------------------------------------
    # # Abslute value on left, flat on right.
    # f = lambda x: np.where(x <= 0, abs(x), 0)
    # x = np.linspace(-1.1, 1, 10)
    # y = f(x)

    # # --------------------------------------------------------------------
    # # Sin function.
    # np.random.seed(1)
    # x = np.random.random(size=(30,)) * 6*np.pi
    # x.sort()
    # y = np.sin(x)


    # Watson test
    # random test
    # 

    # - >=3 changes in monotonicity
    # - strictly monotone data
    # - equally spaced x values, unequally spaced x values
    # - isolated, two consecutive, three consecutive
    # - straight lines (>=3 points)
    # - large tangents
    # - points near extreme points


    # We should have some data with several (say >=3) changes in
    # monotonicity, as well as data that is strictly monotone.  We
    # should have examples with equally spaced x values and with very
    # unequally spaced x values.  We should have every possible
    # pattern of extrema (isolated, two consecutive, three
    # consecutive), some data with straight line segments (>=3 points
    # in a row), and some data with very large tangents (e.g., f(x) =
    # 1+(1/(x-1)) with x_i from 0 to .98.  And any other "edge" case
    # you can think of. 


    # Plot the new test case and the approximation performance.

    if False:
        # Generate a single example plot (with oversied markers).
        p = example_points(x,y, ms=10)
        p.plot(file_name="pmqsi_example.html", show=False)

    if True:
        # Generate individual application plots, with normal sized markers.
        # visualize_one(x, y, cubic=True, show=True)
        visualize_one(x, y, cubic=False, show=True, verbose=True)


# Find places of "agreement in change" by looking at the harmonic mean
# of quadratic slopes coming from the left and right.
# 
# Start with the location that has the greatest difference between
# current slope and target slope.
# 
# Fix that point by doing a binary search between current Df DDf, and
# the target Df DDf, landing on the closest point to the target that
# is still monotone.
# 
# Re-evaluate the estimations for the neighbors of the fixed point.
# 
# Fix the next point in the list until all points have been fixed once.
# 

