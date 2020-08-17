
# Custom exception for when there are data errors.
class BadUsage(Exception): pass
class BadData(Exception): pass
class BadValues(Exception): pass
class UnexpectedError(Exception): pass

# Given the left and right slope ratios (alpha and beta in paper),
# compute a new left and right slope ratio that rest in the feasible
# region for a monotonic cubic interpolating polynomial.
# 
#   mode -- 1 (curvy)   Project onto the box inside of "3".
#           2 (default) Project onto circle of radius "3".
#           3 ( ~~ )    Project onto line "y = 3 - x".
#           4 (linear)  Project onto max of "y = 3 - 2x" and "x = 3 - 2y".
# 
#   project -- 1 uses the intercept of a line with the region.
#              2 uses the closest point on the region boundary.
# 
def compute_cubic_feasible(left_ratio, right_ratio, mode, project):
    # First check to make sure we have valid ratios (no division by 0).
    if (max(left_ratio, right_ratio) == 0):
        pass
    # Project the ratios onto the appropriate surface.
    elif (mode == 1):
        # Option 1: Project onto the square of side length 3.
        if project == 1:
            # Determine which line segment it will project onto.
            mult = 3 / max(left_ratio, right_ratio)
            # Perform the projection onto the line segment by
            # shortening the vector to make the max coordinate 3.
            left_ratio = min(left_ratio, mult*left_ratio)
            right_ratio = min(right_ratio, mult*right_ratio)
        elif (project == 2):
            # The projection direction is naively determined by
            # capping the value of the vector coordinates at 3.
            left_ratio = min(left_ratio, 3)
            right_ratio = min(right_ratio, 3)
    elif (mode == 2):
        # Option 2: Project onto the circle with radius 3.
        #           Same behavior regardless of projection method.
        length = (left_ratio**2 + right_ratio**2)**(1/2)
        if (length > 3):
            left_ratio /= (length / 3)
            right_ratio /= (length / 3)
    elif (mode == 3):
        # Option 3: Project onto the line "y = 3 - x".
        if (right_ratio != 0): slope = left_ratio / right_ratio
        else:                  slope = float('inf')
        if (project == 1):
            # Compute the intersection of the line and the ratio vector.
            left_ratio = min(left_ratio, (slope*3) / (slope + 1))
            right_ratio = min(right_ratio, 3 / (slope + 1))
        elif (project == 2):
            # Compute the projection onto the line.
            left_ratio, right_ratio = (
                min(3, left_ratio,  max(0,(left_ratio - right_ratio + 3)/2)),
                min(3, right_ratio, max(0,(right_ratio - left_ratio + 3)/2))
            )
    elif (mode == 4):
        # Option 4: Project onto the greater of "y = 3 - 2x" and "x = 3 - 2y".
        if (right_ratio != 0): slope = left_ratio / right_ratio
        else:                  slope = float('inf')
        if (project == 1):
            # Intersect the ratio vector with the first (steeper) line.
            first_left = min(left_ratio, (3*slope) / (2*slope + 1))
            first_right = min(right_ratio, 3 / (2*slope + 1))
            # Intersect the ratio vector with the second (shallower) line.
            second_left = min(left_ratio, (3*slope) / (slope + 2))
            second_right = min(right_ratio, 3 / (slope + 2))
            # Compute the intersection of the line and the ratio vector.
            left_ratio = min(left_ratio, max(first_left, second_left))
            right_ratio = min(right_ratio, max(first_right, second_right))
        elif (project == 2):
            if   (left_ratio >= right_ratio):
                # Compute the projection onto the first (steeper) line.
                left_ratio, right_ratio = (
                    min(3, left_ratio,  max(0,(3 + 4*left_ratio - 2*right_ratio)/5)),
                    min(3, right_ratio, max(0,(6 - 2*left_ratio +   right_ratio)/5)) )
            elif (right_ratio >  left_ratio):
                # Compute the projection onto the second (shallower) line.
                left_ratio, right_ratio = (
                    min(3, left_ratio,  max(0,(6 +   left_ratio - 2*right_ratio)/5)),
                    min(3, right_ratio, max(0,(3 - 2*left_ratio + 4*right_ratio)/5)) )
            # else:
            #     # Project the point exactly onto the intersection.
            #     left_ratio, right_ratio = 1, 1
    return left_ratio, right_ratio


# Given data values and function values for a monotone function over
# the interval defined by [min(x), max(x)], construct a piecewise
# cubic monotone interpolating spline function of "x" and return it.
# 
#   mode -- 1 (curvy)   Project onto the box inside of "3".
#           2 (default) Project onto circle of radius "3".
#           3 ( ~~ )    Project onto line "y = 3 - x".
#           4 (linear)  Project onto max of "y = 3 - 2x" and "x = 3 - 2y".
# 
#   ends -- 0 locks derivatives at ends to 0.
#           1 locks derivatives at ends to be secant slope to neighbor.
#           2 locks derivatives at ends to be quadratic slope through 2 neighbors.
#           ~ manually specify a 2-tuple that declares the endpoint derivatives.
# 
#   mids -- 0 locks derivatives at midpoints to 0.
#           1 locks derivatives at midpoints to be secant slope through neighbors.
#           2 locks derivatives at midpoints to be quadratic slope through neighbors.
#           ~ manually specify an (n-2)-tuple that declares the midpoint derivatives.
# 
#   project -- 1 uses the intercept of a line with the region.
#              2 uses the closest point on the region boundary.
# 
def cubic_interp(x, y, mode=2, ends=0, mids=2, project=1):
    # ----------------------------------------------------------------
    # STEP 0: Check for proper usage, organize data.
    # 
    # Verify the lengths of the data values and function values are the same.
    if (len(x) != len(y)):
        raise(BadUsage(f"Received '{len(x)}' points and '{len(y)}' function values, should be same number."))
    # Make sure the data is sorted correctly (by data values).
    sort_order = sorted(range(len(x)), key=lambda i: x[i])
    if (sort_order != list(range(len(x)))):
        x = [x[i] for i in sort_order]
        y = [y[i] for i in sort_order]
    # Check that data is non-repeating.
    if (sorted(set(x)) != x):
        raise(BadData("All data values must be unique, duplicated data values are not permitted."))
    # Check that the data is monotone.
    non_decreasing = all(y[i] <= y[i+1] for i in range(len(y)-1))
    non_increasing = all(y[i+1] <= y[i] for i in range(len(y)-1))
    if not (non_increasing or non_decreasing):
        raise(BadValues("Received function values that are not monotone."))
    # ----------------------------------------------------------------
    # STEP 1: Initialize all derivatives to reasonable values.
    #  
    # ends:
    #   (zero)   endpoints to zero.
    #   (lin)    endpoints to secant slope through endpoint neighbor.
    #   (quad)   endpoints to capped quadratic interpolant slope.
    #   (manual) a 2-tuple provides locked-in values for derivatives.
    # 
    # mids:
    #   (zero)   all slopes are locked into the value 0.
    #   (lin)    secant slope between left and right neighbor.
    #   (quad)   slope of quadratic interpolant over three point window.
    #   (manual) an (n-2)-tuple provides locked-in values for derivatives.
    #
    from polynomial import fill_derivative
    deriv = fill_derivative(x, y, ends, mids)
    # ----------------------------------------------------------------
    # STEP 2: Fix the derivatives that violate monotonicity.
    # 
    #         We have to be careful though, our procedure must not
    #         break monotonicity on a previous interval while
    #         fixing the current interval.
    # 
    #         We will analyze the ratios of assigned derivatives
    #         to the slope of the secant line between data points.
    #         Respectively the coefficients "alpha" and "beta"
    #         from the original 1980 paper.
    # 
    #  1 (curvy)   Project onto the box inside of "3".
    #  2 (default) Project onto circle of radius "3".
    #  3 ( ~~ )    Project onto line "y = 3 - x".
    #  4 (linear)  Project onto max of "y = 3 - 2x" and "x = 3 - 2y".
    # 
    #         In the paper, the projections are done by
    #         intersecting a line through the origin and the
    #         the "alpha, beta" point defined by current
    #         derivatives.
    # 
    # Looping through derivatives and adjust them to enforce monotonicity.
    for i in range(len(x)-1):
        # Compute the secant slope, the left slope ratio and the
        # right slope ratio for this interval of the function.
        secant_slope = (y[i+1] - y[i]) / (x[i+1] - x[i])
        left_ratio = deriv[i] / secant_slope
        right_ratio = deriv[i+1] / secant_slope
        # Compute the (projected) left and right ratios that ensure
        # the cubic function is monotonic over the interval.
        left_ratio, right_ratio = compute_cubic_feasible(
            left_ratio, right_ratio, mode, project)
        # Set the derivative values based on (projected) monotone slope ratios.
        deriv[i]   = left_ratio * secant_slope
        deriv[i+1] = right_ratio * secant_slope
    # ----------------------------------------------------------------
    # STEP 3: Define an piecewise polynomial interpolating function.
    from polynomial import Spline
    knots = x
    values = [[f,df] for (f,df) in zip(y, deriv)]
    return Spline(knots, values)

# Testing code.
if __name__ == "__main__":
    TEST_FIT = True
    SMALL_FIT_DEMO = False
    TEST_PROJECTION = False
    SHOW_FEASIBLE_REGION = False
    BUILDING_APPROX_DEMO = False
    SMALL_PROJECTION_DEMO = False

    if SHOW_FEASIBLE_REGION:
        # Compute the circle boundary points (for region 2).
        circ_func = lambda x: (9 - x**2)**(1/2)
        x = [(i/100) for i in range(300+1)] + [0]
        y = [circ_func(i) for i in x]

        from util.plot import Plot
        p = Plot("Regions of monotonicity", "alpha", "beta", font_family="times")

        # Add region 1 (the box at 3).
        p.add("Region 1", [0,3,3]+x[::-1][1:], [3,3,0]+y[::-1][1:], mode="lines", fill="toself")
        # Add region 2 (the circle at 3).
        p.add("Region 2", x, y, mode="lines", fill="toself")
        # Add region 3 (inside y = 3 - x).
        p.add("Region 3", [0,3,1,0], [3,0,1,3], mode="lines", fill="toself")
        # Add region 4 (inside min{y = 3 - 2x, x = 3 - 2y}).
        p.add("Region 4", [0,0,1,3,0], [0,3,1,0,0], mode="lines", fill="toself")
        # Define knowing-to-work width and height, compute ratio.
        width, height= 880, 700
        ratio = width / height
        # Redefine the width an height
        width = 500
        height = 400
        p.show(x_range=[-1,4], y_range=[-1,4], width=width,
               height=height, file_name="feasible_region.html")
        

    if SMALL_PROJECTION_DEMO:
        from util.plot import Plot
        # Creat an animated plot.
        p = Plot("Example Projections onto the Region of Monotonicity",
                 "alpha", "beta", font_family="times")

        # Compute the circle boundary points (for region 2).
        circ_func = lambda x: (9 - x**2)**(1/2)
        x = [(i/100) for i in range(300+1)] + [0]
        y = [circ_func(i) for i in x]

        # Add points and their projections.
        def add_point(x, y, color=None, mode=2, project=1):
            show = False
            if (color == None):
                color = p.color_num + 1
                p.color_num += 1
                show = True
            new_y, new_x = compute_cubic_feasible(y, x, mode=mode, project=project)
            p.add(f"Point {color-3} to Region {mode}", [x,new_x], [y,new_y],
                  mode="markers+lines", dash="dot", color=p.color(4),
                  line_color=p.color(mode-1, alpha=.6))
            return color

        # Add region 1 (the box at 3).
        p.add("Region 1", [0,3,3]+x[::-1][1:], [3,3,0]+y[::-1][1:], mode="lines", fill="toself")
        # Add region 2 (the circle at 3).
        p.add("Region 2", x, y, mode="lines", fill="toself")
        # Add region 3 (inside y = 3 - x).
        p.add("Region 3", [0,3,1,0], [3,0,1,3], mode="lines", fill="toself")
        # Add region 4 (inside min{y = 3 - 2x, x = 3 - 2y}).
        p.add("Region 4", [0,0,1,3,0], [0,3,1,0,0], mode="lines", fill="toself")

        add_point(3,4.5, mode=1)
        add_point(4,4, mode=2)
        add_point(4,1, mode=3)
        add_point(1.5,3.5, mode=4)
        p.show(x_range=[-1,5], y_range=[-1,5], width=500,
               height=400, file_name="demo_projection.html")


    if SMALL_FIT_DEMO:
        from util.plot import Plot
        SEED = 1
        NODES = 4

        # Generate random data to test the monotonic fit function.
        import numpy as np
        np.random.seed(SEED)
        nodes = NODES + 2
        x = np.linspace(0, 1, nodes)
        y = sorted(np.random.normal(size=(nodes,)))
        # Convert these arrays to lists.
        x, y = list(x), list(y)

        # Generate a plot of the function.
        width = max(x) - min(x)
        padding = width * .1
        bounds = (min(x)-padding, max(x)+padding)
        # Reserve the "0" color for the points (to be added last).
        p = Plot("Cubic Montone Interpolating Spline", font_family="times")
        p.color_num += 1
        # Settings for the  custom fits.
        project = 1
        ends = 1
        mids = 2
        mode = 2
        fit = cubic_interp(x, y, mode=mode, mids=mids,
                           project=project, ends=ends)
        p.add_func(f"Fit: mode {mode}, project {project}, ends {ends}, mids {mids}", 
                   fit, bounds, color=p.color(mode+2))
        # Add the points last (so they are on top).
        p.add("Points", x, y, color=p.color(0))
        p.show(file_name="demo_fit.html", show_legend=False,
               width=500, height=400)

        p = Plot("Derivative of Cubic Montone Spline", font_family="times")
        p.color_num += 1
        fit = cubic_interp(x, y, mode=mode, mids=mids,
                           project=project, ends=ends)
        for (x1, y1) in zip(x,y):
            p.add("Point", [x1, x1], [0,9], mode="lines",
                  color=p.color(0, alpha=.3))

        p.add_func(f"DERIVATIVE Fit: mode {mode}, project {project}, ends {ends}, mids {mids}", 
                   fit.derivative(1), bounds, color=p.color(mode+2))

        p.show(file_name="demo_fit_deriv.html", show_legend=False,
               width=500, height=400)


    if TEST_PROJECTION:
        mode = 2
        project = 1

        # Compute the circle boundary points (for region 2).
        circ_func = lambda x: (9 - x**2)**(1/2)
        x = [(i/100) for i in range(300+1)] + [0]
        y = [circ_func(i) for i in x]

        for m in (1,2,3,4):
            for pr in (1, 2):
                from util.plot import Plot
                # Creat an animated plot.
                p = Plot("Projecting onto the region of monotonicity", "d[i+1] / s", "d[i] / s")

                # Add points and their projections.
                def add_point(x, y, color=None, mode=2, project=2):
                    show = False
                    if (color == None):
                        color = p.color_num + 1
                        p.color_num += 1
                        show = True
                     new_y, new_x = compute_cubic_feasible(y, x, mode=mode, project=project)
                    p.add(f"Point {color-4} to Region {mode}, {project}", [x,new_x], [y,new_y],
                          mode="markers+lines", dash="dot", color=p.color(4),
                          line_color=p.color(mode-1, alpha=.3))
                    return color

                # Add region 1 (the box at 3).
                p.add("Region 1", [0,3,3]+x[::-1][1:], [3,3,0]+y[::-1][1:], mode="lines", fill="toself")
                # Add region 2 (the circle at 3).
                p.add("Region 2", x, y, mode="lines", fill="toself")
                # Add region 3 (inside y = 3 - x).
                p.add("Region 3", [0,3,1,0], [3,0,1,3], mode="lines", fill="toself")
                # Add region 4 (inside min{y = 3 - 2x, x = 3 - 2y}).
                p.add("Region 4", [0,0,1,3,0], [0,3,1,0,0], mode="lines", fill="toself")

                c = add_point(.3,.3, color=None,  mode=m, project=pr)
                c = add_point(1,1, color=None,  mode=m, project=pr)
                c = add_point(2,2, color=None,  mode=m, project=pr)
                c = add_point(3,3, color=None,  mode=m, project=pr)
                c = add_point(4,4, color=None,  mode=m, project=pr)

                c = add_point(2.1,2, color=None,  mode=m, project=pr)
                c = add_point(3.1,3, color=None,  mode=m, project=pr)
                c = add_point(4.1,4, color=None,  mode=m, project=pr)

                c = add_point(2,4, color=None,  mode=m, project=pr)
                c = add_point(2,5, color=None,  mode=m, project=pr)
                c = add_point(2,8, color=None,  mode=m, project=pr)
                c = add_point(4,8, color=None,  mode=m, project=pr)
                c = add_point(4,8, color=None,  mode=m, project=pr)
                c = add_point(4,4.3, color=None,  mode=m, project=pr)

                c = add_point(4,2, color=None,  mode=m, project=pr)
                c = add_point(5,2, color=None,  mode=m, project=pr)
                c = add_point(8,2, color=None,  mode=m, project=pr)
                c = add_point(8,4, color=None,  mode=m, project=pr)
                c = add_point(4.3,4, color=None,  mode=m, project=pr)

                append = not ((m == 1) and (pr == 1))
                p.show(x_range=[-1,9], y_range=[-1,9], width=880,
                       height=700, append=append, file_name="demo_feasible.html")


    if TEST_FIT:
        from util.plot import Plot
        SEED = 1
        NODES = 10
        QUADS = False
        FUNCS = True
        DERIVS = True

        # Generate random data to test the monotonic fit function.
        import numpy as np
        np.random.seed(SEED)
        nodes = NODES + 2
        x = np.linspace(0, 1, nodes)
        # y = [0.]
        # for i in range(1,nodes-1):
        #     new_val = y[-1] + np.random.random()**(nodes-i-1) * (1. - y[-1])
        #     while (new_val == y[-1]):
        #         new_val = y[-1] + np.random.random()**(nodes-i-1) * (1. - y[-1])
        #     y.append( new_val )
        # y.append(1.)
        # y = np.array(y)
        y = sorted(np.random.normal(size=(nodes,)))
        # Convert these arrays to lists.
        x, y = list(x), list(y)
        print("x: ",x)
        print("y: ",y)

        # Generate a plot of the function.
        p = Plot("Cubic montone fits")
        width = max(x) - min(x)
        padding = width * .1
        bounds = (min(x)-padding, max(x)+padding)
        # Reserve the "0" color for the points (to be added last).
        p.color_num += 1
        if QUADS:
            from polynomial import solve_quadratic
            # Add the quadratic fits about pairs of three points.
            for i in range(1, len(x)-1):
                a,b,c = solve_quadratic(x[i-1:i+2], y[i-1:i+2])
                f = lambda x: a*x**2 + b*x + c
                p.add_func(f"{i} quad", f, bounds, color=p.color((0,0,0,.1)), dash="dash")
        # Add the "correct" PCHIP fit.
        from scipy.interpolate import PchipInterpolator
        pchip_fit = PchipInterpolator(x, y)
        # Settings for the  custom fits.
        project = 1
        ends = 0
        mids = 1
        # ------------------------------------------------------------
        #    Add all of the interpolating functions.
        if FUNCS:
            p.add_func("PCHIP", pchip_fit, bounds)
            # Add the different fits from my function.
            for mode in (1,2,3,4):
                fit = cubic_interp(x, y, mode=mode, mids=mids,
                                   project=project, ends=ends)
                p.add_func(f"Fit: mode {mode}, project {project}, ends {ends}, mids {mids}", 
                           fit, bounds, color=p.color(mode+2))
        # ------------------------------------------------------------
        #    Add all of the derivatives of the interplating functions.
        if DERIVS:
            p.add_func("DERIVATIVE PCHIP", pchip_fit.derivative(1), bounds)
            # Add the different fits from my function.
            for mode in (1,2,3,4):
                fit = cubic_interp(x, y, mode=mode, mids=mids,
                                   project=project, ends=ends)
                p.add_func(f"DERIVATIVE Fit: mode {mode}, project {project}, ends {ends}, mids {mids}", 
                           fit.derivative(1), bounds, color=p.color(mode+2))

        # Add the points last (so they are on top).
        p.add("Points", x, y, color=p.color(0))
        p.show(file_name="demo_fit_many.html")


    if BUILDING_APPROX_DEMO:
        from util.plot import Plot
        from util.stats import cdf_fit, cdf_points
        import numpy as np
        # Generate random data.
        np.random.seed(4)
        data = np.random.normal(size=(8,))
        data = sorted(data + min(data))
        print([round(v,2) for v in data])
        x, y = cdf_points(data)
        x[0] -= .1

        # # Make a flat EDF plot.
        # p = Plot()
        # p.add("EDF Points", x, y)
        # f = cdf_fit((x,y), fit=None)
        # p.add_func("Flat fit", f, (min(data)-.2, max(data)+.2))
        # p.show(file_name="flat-fit.html")

        # # Make a linear fit.
        # p = Plot()
        # p.add("EDF Points", x, y)
        # f = cdf_fit((x,y), fit="linear")
        # p.add_func("Linear fit", f, (min(data)-.2, max(data)+.2))
        # p.show(file_name="linear-fit.html")

        # # Make a bad cubic fit.
        # p = Plot()
        # p.add("EDF Points", x, y)
        # f = cdf_fit((x,y), fit="cubic")
        # p.add_func("Cubic fit", f, (min(data)-.2, max(data)+.2))
        # p.show(file_name="cubic-fit-bad.html")

        # Make a good cubic fit.
        p = Plot()
        p.add("EDF Points", x, y)
        f = cubic_interp(list(x),list(y), mode=4, mids=2)
        p.add_func("Cubic fit", f, (min(data)-.2, max(data)+.2))
        p.show(file_name="cubic-fit-good.html")

        # f = cdf_fit((x,y), fit="linear")
        # p.add_func("Linear fit", f, (min(data)-.2, max(data)+.2))
        # f = cdf_fit((x,y), fit="cubic")
        # p.add_func("Cubic fit", f, (min(data)-.2, max(data)+.2))

