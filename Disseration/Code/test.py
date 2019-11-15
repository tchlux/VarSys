import numpy as np
import fmodpy
from util.plot import Plot

splines = fmodpy.fimport("splines.f08", verbose=True,
                         autocompile_extra_files=True,
                         module_link_args=["-lblas", "-llapack"])


TEST_BANDED_MATRIX = False
TEST_B_SPLINE = True
TEST_SPLINE = False
TEST_SPLINE_L2 = False


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
    knots = np.array([0.,1,2,3,4])
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
