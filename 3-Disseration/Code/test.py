import numpy as np
import fmodpy
from util.plot import Plot

splines = fmodpy.fimport("splines.f08", verbose=True,
                         autocompile_extra_files=True,
                         module_link_args=["-lblas", "-llapack"])


TEST_BANDED_MATRIX = False
TEST_B_SPLINE = True
TEST_SPLINE = False
TEST_SPLINE_ADD = False
TEST_SPLINE_L2 = False
TEST_MONOTONE = False

if TEST_MONOTONE:
    # Return a function that interpolates the given values at the given knots.
    def fit_spline(x, y):
        # Make sure "y" is two dimensional.
        try: len(y[0])
        except: y = [[v] for v in y]
        # Fit the spline and get the knots / coefficients for the
        # B-spline representation of the data.
        sk, sc = splines.fit_spline(np.array(x, order='F', dtype=float), 
                                    np.array(y, order='F', dtype=float))
        # Define an interpolating function with some stored attributes.
        def f(x, d=0):
            try:    len(x)
            except: x = np.array([x], order="f", dtype=float)
            return splines.eval_spline(sk, sc, x.copy(), d=d)
        f.x = x.copy()
        f.y = y.copy()
        f.knots = sk
        f.coefs = sc
        return f


    SEED = 0
    NODES = 30
    SUBSET = slice(None)

    # Generate random data to test the monotonic fit function.
    import numpy as np
    np.random.seed(SEED)
    nodes = NODES + 2
    x = np.linspace(0, 1, nodes)
    y = sorted(np.random.normal(size=(nodes,)))
    x -= min(x); x /= max(x)
    y -= min(y); y /= max(y); y -= max(y)
    # Convert these arrays to lists.
    x, y = list(x)[SUBSET], list(y)[SUBSET]
    # Convert these arrays to exact arithmetic.
    from util.math import Fraction
    x = list(map(Fraction, x))
    y = list(map(Fraction, y))
    interval = [float(min(x)), float(max(x))]
    from polynomial import fit
    f = fit(x, y, continuity=1)
    x = np.array([float(v) for v in f.knots])
    y = np.array([[float(v) for v in l] for l in f.values])
    b = fit_spline(x, y)
    # Make a plot.
    p = Plot()
    p.add("points", x, y[:,0])
    p.add_func("f", f, [min(x), max(x)])
    p.add_func("b", b, [min(x), max(x)])

    # max_val = -float('inf')
    # for i in range(0, len(b.coefs)):
    #     max_val = max(b.coefs[i], max_val)
    #     b.coefs[i] = max_val

    min_val = float('inf')
    for i in range(len(b.coefs)-1,-1,-1):
        min_val = min(min_val, b.coefs[i])
        b.coefs[i] = min_val

    p.add_func("b monotone", b, [min(x), max(x)])

    # e_coefs = b.coefs.copy() * 0
    # for _ in range(1):
    #     # Get the errors, but clip all derivatives to be positive.
    #     error = y.copy()
    #     for i in range(y.shape[1]):
    #         error[:,i] = y[:,i] - b(x, d=i)
    #     np.clip(error[:,1], 0, None, out=error[:,1])
    #     # Fit a B-spline fit to the error terms.
    #     e = fit_spline(x, error)
    #     # Make the B-spline fit monotone.
    #     min_val = float('inf')
    #     for i in range(len(e.coefs)-1,-1,-1):
    #         min_val = min(min_val, e.coefs[i])
    #         e.coefs[i] = min_val
    #     # Add the B-spline correction to the monotone version.
    #     b.coefs += e.coefs
    #     # Sum the total of the corrections that have been made.
    #     e_coefs += e.coefs

    # # Make the "e" function the sum of all corrections.
    # e.coefs[:] = e_coefs[:]
    # p.add_func("correction", e, [min(x), max(x)])

    # # Show off the correct "b" function.
    # p.add_func("b+corrections", b, [min(x), max(x)])

    evals = y.copy()
    for i in range(y.shape[1]):
        evals[:,i] = b(x, d=i)
    print()
    print("y:")
    print(y)
    print()
    print("b:")
    print(evals)
    print()
    print('errors:')
    print(y - evals)
    print()

    p.show(file_name="monotone-b-spline.html")




if TEST_SPLINE_ADD:
    # Return a function that interpolates the given values at the given knots.
    def fit_spline(x, y):
        # Make sure "y" is two dimensional.
        try: len(y[0])
        except: y = [[v] for v in y]
        # Fit the spline and get the knots / coefficients for the
        # B-spline representation of the data.
        sk, sc = splines.fit_spline(np.array(x, order='F', dtype=float), 
                                    np.array(y, order='F', dtype=float))
        # Define an interpolating function with some stored attributes.
        def f(x, d=0):
            try:    len(x)
            except: x = np.array([x], order="f", dtype=float)
            return splines.eval_spline(sk, sc, x.copy(), d=d)
        f.x = x.copy()
        f.y = y.copy()
        f.knots = sk
        f.coefs = sc
        return f



    # K1 size:  7  K2 size: 16
    # K1:   1.00  1.00  1.00  1.50  2.00  3.00  3.00
    # K2:   1.25  1.25  1.25  1.25  1.50  1.50  1.75  1.75  2.00  2.00  2.50  2.50  5.00  5.00  5.00  5.00
    #
    # 1 I1:   1   I2:   1  1.00  1.25  0.00
    # 2 I1:   4   I2:   1  1.50  1.25  0.00
    # 1 I1:   4   I2:   7  1.50  1.75  0.00
    # 2 I1:   5   I2:   7  2.00  1.75  0.00
    # 1 I1:   5   I2:  11  2.00  2.50  0.00
    # 2 I1:   6   I2:  11  3.00  2.50  0.00
    # 1 I1:   6   I2:  13  3.00  5.00  0.00


    knots1 = [1, 1, 1.5, 2, 3]
    values1 = [(i,) for i in range(len(knots1))]
    f1 = fit_spline(knots1, values1)

    knots2 = [1.25, 1.5, 1.75, 2, 2.5, 5]
    values2 = [(i,0) for i in range(len(knots2))]
    f2 = fit_spline(knots2, values2)

    p = Plot()
    p.add_func("F1", f1, [0, 5])
    p.add_func("F2", f2, [0, 5])
    p.add_func("F1d", lambda x: f1(x,d=1), [0, 5])
    p.show()
    exit()
    print("knots1: ",len(knots1), knots1)
    print("knots2: ",len(knots2), knots2)
    ans = sorted(set(knots1).union(knots2))
    print("answer: ", len(ans), ans)
    sk, sc, nk, nc = splines.add_splines(f1.knots, f1.coefs, f2.knots, f2.coefs)


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
    # knots = np.array([0.,1,1,1])
    knots = np.array([0.,1,2,3])
    x = np.linspace(min(knots)-.1,max(knots)+.1,101)
    y = splines.eval_bspline(knots, x.copy(), d=0)

    from scipy.integrate import quad as integrate
    # Evaluate the numerical derivative (or integral) of a B-spline.
    def truth(x, d=0, f=None):
        if (d == 0): return splines.eval_bspline(knots,np.array([x]))[0]
        if (d > 0):
            return (truth(x+2**(-8), d-1) - truth(x-2**(-8), d-1)) / 2**(-7)
        if (d < 0):
            # Initialize the function to be integrated, if not given.
            if (f is None): f = truth
            # If this is a raw integral, just compute it directly.
            if (d == -1): return integrate(f, knots[0], x)[0]
            else:
                # Otherwise integrate the given function (and recurse).
                return truth(x, d=d+1, f=lambda x: integrate(f, knots[0], x)[0])

    print()
    for d in range(len(knots)-1):
        at_knots = splines.eval_bspline(knots,knots.copy(), d=d)
        at_true  = np.array([truth(v,d) for v in knots])
        print(f"D{d} at knots: ", at_knots)

    for d in range(1,3):
        at_knots = splines.eval_bspline(knots,knots.copy(), d=-d)
        at_true  = np.array([truth(v,-d) for v in knots])
        print()
        print(f"I{d} at knots: ", at_knots)
        print(f"I{d} truth:    ", at_true)
        print(f"I{d}           ", at_true / at_knots)
    print()

    # Make a pretty picture of the B-spline and its derivatives.
    styles = ["dot", "dash", "dashdot"]
    # Add the B-spline function itself.
    continuity = "" if (len(knots) <= 2) else f"C{len(knots)-3} "
    p = Plot(continuity+f"B-Spline, piecewise order {len(knots)-1} polynomial")
    p.add("B-spline", x, y, mode="markers+lines", color=p.color(1))
    # Add all interesting derivatives.
    for d in range(1, len(knots)-1):
        color = p.color(d-1) if d == 1 else p.color(d)
        dash = styles[(d-1)%len(styles)]
        dy = splines.eval_bspline(knots,x.copy(), d=d)
        name = f"D{d}"
        p.add(name, x, dy, mode="lines", color=color, dash=dash, group=name)
        color = p.color(d-1,alpha=.3) if d == 1 else p.color(d,alpha=.3)
        p.add_func(f"True {name}", lambda x: truth(x, d=d), [knots[0],knots[-1]], color=color, group=name)

    # Add all interesting derivatives.
    for d in range(1,3):
        i = d + len(knots)-2
        color = p.color(i)
        dash = styles[(i-1)%len(styles)]
        dy = splines.eval_bspline(knots,x.copy(), d=-d)
        name = f"I{d}"
        p.add(name, x, dy, mode="lines", color=color, dash=dash, group=name)
        color = p.color(i, alpha=.3)
        p.add_func(f"True {name}", lambda x: truth(x, d=-d), [knots[0],knots[-1]], color=color, group=name,plot_points=10)

    # Add the knots themselves.
    p.add("Knot Locations", knots, [0]*len(knots), color=p.color(0))

    # Generate the visual.
    p.show(file_name=f"deriv_test-{len(knots)-1}.html", y_range=[-2,5.5])


# ==================================
#      Test the spline-fit code     
# ==================================

if TEST_SPLINE:
    VISUALIZE_TEST = True

    if VISUALIZE_TEST:
        knots = np.array([0.,1,2])
        # values = np.array([[1.,0,2,1], [.5,0,1,2], [0,0,0,3]])
        values = np.array([[1.,0], [0.,0], [0.,0]])
        continuity = values.shape[1]-1
    else:
        num_knots = 100
        continuity = 4
        # Generate some random knots and values.
        knots = np.random.random(size=(num_knots))
        knots.sort()
        values = np.random.random(size=(num_knots,continuity))

    # Get the spline knots and spline coefficients.
    values = np.asfortranarray(values)
    sk, sc = splines.fit_spline(knots, values)
    sc *= -1
    for d in range(continuity):
        y = splines.eval_spline(sk, sc, knots.copy(), d=d)
        error = abs(y - values[:,d])
        print(f"f^({d}) max error: {max(error):.3e}")

    if VISUALIZE_TEST:
        print()
        print("sk: ",sk)
        print("sc: ",sc)
        print()
        padding = 0 #.1
        x = np.linspace(min(knots)-padding,max(knots)+padding,(len(knots)-1)*100+1)
        y = splines.eval_spline(sk, sc, x.copy(), d=0)

        # Make a pretty picture of the B-spline and its derivatives.
        p = Plot("Hermite Interpolant")
        p.add("Spline", x, y, mode="lines", color=p.color(1), group="s")
        p.add("Knots", knots, values[:,0], color=p.color(1), group="s")
        # Add all interesting derivatives.
        styles = ["dot", "dash", "dashdot"]
        for d in range(1, values.shape[1]):
            color = p.color(d-1) if d == 1 else p.color(d)
            dash = styles[(d-1)%len(styles)]
            dy = splines.eval_spline(sk, sc, x.copy(), d=d)
            p.add(f"D{d}", x, dy, mode="lines", color=color, dash=dash,group=d)
            p.add(f"k{d}", knots, values[:,d], color=color, group=d)
        y = splines.eval_spline(sk, sc, x.copy(), d=-1)
        p.add("Integral", x, y, mode="lines", color=p.color(0), group="i")

        p.show(file_name=f"spline_test-N{len(values)}-C{len(values[0])}.html")
 


if TEST_SPLINE_L2:

    # Add a plot of a spline object.
    def plot_spline(name, f, p=Plot()):
        knots = f.x
        values = f.y
        x = np.linspace(min(knots),max(knots),(len(knots)-1)*100+1)
        y = f(x)
        p.add(f"{name} spline", x, y, mode="lines", color=p.color(1), group=f"{name}s")
        p.add(f"{name} knots", knots, values[:,0], color=p.color(1), group=f"{name}s")
        # Add all interesting derivatives.
        styles = ["dot", "dash", "dashdot"]
        for d in range(1, values.shape[1]+1):
            color = p.color(d-1) if d == 1 else p.color(d)
            dash = styles[(d-1)%len(styles)]
            dy = f(x, d=d)
            p.add(f"{name} D{d}", x, dy, mode="lines", color=color, dash=dash,group=f"{name}{d}")
            if (d < values.shape[1]):
                p.add(f"{name} D{d} knots", knots, values[:,d], color=color, group=f"{name}{d}")
        # Return the plot object.
        return p

    # Compute the numerical approximation to the derivative of 'f' at 'x'.
    def numerical_derivative(f, x, d=1, width=2**(-13)):
        if (d == 1):
            return (f(x + width) - f(x)) / width
        else:
            return (numerical_derivative(f, x+width, d=d-1, width=width) -
                    numerical_derivative(f, x, d=d-1, width=width)) / width

    # Return a function that interpolates the given values at the given knots.
    def spline_fit(x, y):
        sk, sc = splines.fit_spline(np.asfortranarray(x), np.asfortranarray(y))
        # Define an interpolating function with some stored attributes.
        def f(x, d=0):
            return splines.eval_spline(sk, sc, x.copy(), d=d)
        f.x = x.copy()
        f.y = y.copy()
        f.knots = sk
        f.coefs = sc
        return f

    # Given some values, compute the values for the squared version of the function.
    def compute_square(knots, values):
        sq_values = []
        # Construct a fit over the knots and values.
        fit = spline_fit(knots, values)
        def fit_squared(x):
            try:    len(x)
            except: x = np.array([x], dtype=float)
            return (splines.eval_spline(fit.knots, fit.coefs, x.copy()))**2
        # Compute the values and derivatives directly.
        for k in range(len(values)):
            new_vals = [
                values[k][0]**2,
            ]
            # Add the first derivative if it applies.
            if len(values[k]) > 1:
                new_vals.append( 2 * (values[k][0] * values[k][1]) )
            # Compute the remaining values with an approximation.
            for d in range(min(len(values[k]),2), len(values[k])*2):
                v = knots[k]
                if (v == knots[-1]): v -= 2**(-8)
                new_vals += [round(numerical_derivative(fit_squared, v, d)[0])]
            # Store the values at this knot.
            sq_values.append(new_vals)
        return np.array(sq_values)

    x = np.array([0, 2], dtype=float)
    # Make a pretty picture of the B-spline and its derivatives.
    p = Plot("Hermite Interpolant")

    y1 = np.array([[1], [2]], dtype=float, order="f")
    # y1 = np.array([[0,1,2], [2,1,0]], dtype=float)
    f = spline_fit(x, y1)
    plot_spline("f", f, p)


    y2 = np.array([[1,0], [2,0]], dtype=float, order="f")
    # y1 = np.array([[0,1,2], [2,1,0]], dtype=float)
    g = spline_fit(x, y2)
    plot_spline("g", g, p)

    print("splines.l2(x,y1,y2): ",splines.l2(x,y1, y2))


    # Compute the:
    #   derivative with respect to coefficients on spline.
    #   derivative with respect to values.


    p.show(file_name=f"spline_difference.html")
    exit()



     # [[  0.  2*  0.   2.  60.]
     # [   4.  2* -2. -54. -12.]]

    print("y1: \n",y1)
    print()

    y3 = compute_square(x, y1)
    # y3[0,2:] = [ 26, -120]
    # y3[1,2:] = [-46, -24]
    print("y3: \n",y3)
    h = spline_fit(x, y3)
    plot_spline("f^2",h,p)

    al = .5

    def f2(x, d=0):
        x = np.array([x], dtype=float)
        return (splines.eval_spline(f.knots, f.coefs, x.copy(), d=d))[0]**2
    p.add_func("f^2", f2, [min(x), max(x)], color=p.color(1, alpha=al))

    d1f2 = lambda x: numerical_derivative(f2, x, d=1, )
    p.add_func("f^2 D1", d1f2, [min(x), max(x)-.001], color=p.color(0, alpha=al))

    d2f2 = lambda x: numerical_derivative(f2, x, d=2)
    p.add_func("f^2 D2", d2f2, [min(x), max(x)-.001], color=p.color(2, alpha=al))

    # d3f2 = lambda x: numerical_derivative(f2, x, d=3)
    # p.add_func("f^2 D3", d3f2, [min(x), max(x)-.001], color=p.color(3, alpha=al))

    # d4f2 = lambda x: numerical_derivative(f2, x, d=4, width=2**(-8))
    # p.add_func("f^2 D4", d4f2, [min(x), max(x)-.02], color=p.color(4, alpha=al))


    p.show(file_name=f"spline_difference.html")

    # Summing and differencing functions is easy.
    # How hard is it to compute a function squared? It will double
    #  the complexity, right? Need to think more.

    # f^0 =    f^0 * f^0
    # f^1 = 2*(f^0 * f^1)
    # f^2 = 2*(f^0 * f^2) +   (f^1 * f^1)
    # f^3 = 2*(f^0 * f^3) + 2*(f^1 * f^2)
    # f^4 = 2*(f^0 * f^4) + 2*(f^1 * f^3) + (f^2 * f^2)
