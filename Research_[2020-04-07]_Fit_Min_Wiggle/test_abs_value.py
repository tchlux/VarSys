import numpy as np
from polynomial import Spline, fit, local_polynomial, local_quadratic, \
    polynomial, Polynomial
from monotone import monotone_quintic_spline
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
def visualize_multiple_methods(x, y, name=["Local quintic"], title="",
                               method=[6], color=[1], dash=[None],
                               **kwargs):
    # Initialize the plot.
    p = Plot(title=title, font_family="times")
    # Add visuals of all methods.
    for (n,m,c,d) in zip(name, method, color, dash):
        if m is not None:
            f = monotone_quintic_spline(x,y, method=m)
        else:
            from scipy.interpolate import PchipInterpolator
            f = PchipInterpolator(x, y)
        p.add_func(n, f, [min(x), max(x)], color=c, dash=d)
    # Add the points styled according to constraints (last to be on top).
    example_points(x,y, ms=6, p=p)
    # Generate the plot and return.
    min_max = [float(min(y)), float(max(y))]
    width = min_max[1] - min_max[0]
    return p.show(y_range=[min_max[0]-width/8, min_max[1]+width/8], **kwargs)



# Visualize the application of one approximation method.
def visualize_one_method(x, y, append=False, title=None, show=False,
                         file_name="pmqsi_example.html", method='c3',
                         local=True, legend=False,  cubic=False, color=1,
                         pure=False, verbose=False, p=None, dash=None):
    if (p is None): p = Plot(font_family="times")

    # Convert "x" and "y" to Fraction types for higher accuracy.
    x = list(map(Fraction, x))
    y = list(map(Fraction, y))

    # Add the cubic interpolant.
    if cubic:
        if dash is None:
            if method is None: file_name, dash = "pchip_example.html", None
            else:                         dash = "dash"
        from scipy.interpolate import PchipInterpolator
        f = PchipInterpolator(x, y)
        p.add_func("Cubic Interpolant", f, [min(x), max(x)], dash=dash, color=2)

    # Add the quintic interpolant.
    if method is not None:
        # Automatically handle the file name and title base on the method.
        file_name = str(method) + "_" + file_name
        m = method
        if (title is not None): title += "constrained"
        if (6 == m):
            if (title is not None): title += " quintic interpolants"
        elif ("r" == m):
            if (title is not None): title += " quadratic regressions"
        elif (4 == m):
            if (title is not None): title += " cubic regressions"
        elif (3 == m):
            if (title is not None): title += " quadratic regressions"
        elif ('facet' in m):
            if (title is not None): title += " local facets"
        else:
            raise(Exception(f"Unknown method '{method}'."))

        # Create the spline.
        local_fits = {}
        f = monotone_quintic_spline(x,y, local_fits=local_fits,
                                    verbose=verbose, method=method)
        p.add_func("Quintic Interpolant", f, [min(x), max(x)], color=color, dash=dash)

        # Add the local interpolants used to generate derivatives.
        if local:
            for i in range(len(x)):
                f, bounds = local_fits[i]
                if (i > 0): bounds[0] = (x[i-1] + x[i]) / 2
                if (i+1 < len(x)):  bounds[1] = (x[i+1] + x[i]) / 2
                p.add_func(f" local approx about {i}", f, bounds, dash="dot", opacity=0.4)

    # Add a pure "fit" that maintains desired continuity.
    if pure:
        p.add_func("pure quintic fit", fit(x, y, continuity=2), [min(x), max(x)], dash="dash")

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

    methods = ['facet', None]
    names = ["facet", "PCHIP"]
    colors = [1, 3]
    dashes = [None, "dot"]
    legend = False
    show = True

    test = -1
    # --------------------------------------------------------------------
    # Watson test.
    y = [0,1,1,1,0,20,19,18,17,0,0,1,0,1]
    x = list(range(len(y)))

    test += 1
    visualize_multiple_methods(x, y, name=names, method=methods,
                               color=colors, dash=dashes,
                               show_legend=legend, show=show,
                               title=f"Test {test}",
                               file_name=f"test_{test}.html")

    # --------------------------------------------------------------------
    # Absolute value function.
    x = np.linspace(-1,1,5)
    y = abs(x)

    test += 1
    visualize_multiple_methods(x, y, name=names, method=methods,
                               color=colors, dash=dashes,
                               show_legend=legend, show=show,
                               title=f"Test {test}",
                               file_name=f"test_{test}.html")

    # --------------------------------------------------------------------
    # Absolute value function.
    x = np.linspace(-1,1,9)
    y = abs(x)

    test += 1
    visualize_multiple_methods(x, y, name=names, method=methods,
                               color=colors, dash=dashes,
                               show_legend=legend, show=show,
                               title=f"Test {test}",
                               file_name=f"test_{test}.html")

    # --------------------------------------------------------------------
    # Abslute value on left, flat on right.
    f = lambda x: np.where(x <= 0, abs(x), 0)
    x = np.linspace(-1.1, 1, 10)
    y = f(x)

    test += 1
    visualize_multiple_methods(x, y, name=names, method=methods,
                               color=colors, dash=dashes,
                               show_legend=legend, show=show,
                               title=f"Test {test}",
                               file_name=f"test_{test}.html")

    # --------------------------------------------------------------------
    # Abslute value on left, flat on right.
    f = lambda x: np.where(x <= 0, abs(x), 0)
    x = np.linspace(-1.1, 1, 20)
    y = f(x)

    test += 1
    visualize_multiple_methods(x, y, name=names, method=methods,
                               color=colors, dash=dashes,
                               show_legend=legend, show=show,
                               title=f"Test {test}",
                               file_name=f"test_{test}.html")

    # --------------------------------------------------------------------
    # Linear function with a small plateau.
    small = 2**(-26)
    x = [-small, 0, 1, 1.9, 2, 2.1, 3, 4, 4+small]
    y = [-small, 0, 1, 1.99, 2, 2.01, 3, 4, 4+small]
    x = x[1:-1]
    y = y[1:-1]

    test += 1
    visualize_multiple_methods(x, y, name=names, method=methods,
                               color=colors, dash=dashes,
                               show_legend=legend, show=show,
                               title=f"Test {test}",
                               file_name=f"test_{test}.html")

    # --------------------------------------------------------------------
    # Linear with one 0 at end.
    x = [0, 1, 2, 6, 7]
    y = [0, 1, 2, 0, -1/2]

    test += 1
    visualize_multiple_methods(x, y, name=names, method=methods,
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
    visualize_multiple_methods(x, y, name=names, method=methods,
                               color=colors, dash=dashes,
                               show_legend=legend, show=show,
                               title=f"Test {test}",
                               file_name=f"test_{test}.html")

    # --------------------------------------------------------------------
    # # Table top level sequence.
    # x = [0, 1, 2, 3, 4]
    # y = [0, 1, 1, 1, 0]



if False:
    # --------------------------------------------------------------------
    # Generate the first visuals explicitly requested by Dr. Watson.

    # Watson test.
    y = [0,1,1,1,0,20,19,18,0,0,1,0,1]
    x = list(range(len(y)))

    if True:
        # Generate a single example plot (with oversied markers).
        p = example_points(x,y, ms=10)
        p.plot(file_name="pmqsi_example.html", show=False)

    # Generate individual application plots, with normal sized markers.
    visualize_one_method(x, y, method=None, cubic=True, show=False)
    visualize_one_method(x, y, method=3, show=False)
    visualize_one_method(x, y, method='facet', show=False)
    visualize_one_method(x, y, method="r", show=False)
