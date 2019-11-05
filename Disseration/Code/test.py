import numpy as np
import fmodpy
from util.plot import Plot

splines = fmodpy.fimport("f_monotone.f08", verbose=True,
                         autocompile_extra_files=True,
                         module_link_args=["-lblas", "-llapack"])


TEST_BANDED_MATRIX = False
TEST_B_SPLINE = False
TEST_SPLINE = True


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
    knots = np.array([0.,1,2,3,4])
    # knots = np.array([0.,1,1,1,1])
    # knots = np.array([0.,0,0,0,1])
    x = np.linspace(min(knots)-.1,max(knots)+.1,101)
    y = splines.eval_bspline(knots,x, d=0)

    print()
    for d in range(len(knots)-1):
        at_knots = splines.eval_bspline(knots,knots, d=d)
        print(f"D{d} at knots: ",at_knots)
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
        dy = splines.eval_bspline(knots,x, d=d)
        p.add(f"D{d}", x, dy, mode="lines", color=color, dash=dash)

    p.show(file_name=f"deriv_test-{len(knots)-1}.html")


# ==================================
#      Test the spline-fit code     
# ==================================

if TEST_SPLINE:
    VISUALIZE_TEST = True

    # knots = np.array([0.,1,2])
    # values = np.array([[1.,0,2,1], [.5,0,1,2], [0,0,0,3]])

    num_knots = 10
    continuity = 5
    # Generate some random knots and values.
    knots = np.random.random(size=(num_knots))
    knots.sort()
    values = np.random.random(size=(num_knots,continuity))
    # Get the spline knots and spline coefficients.
    values = np.asfortranarray(values)
    sk, sv = splines.fit_spline(knots, values)
    for d in range(continuity):
        y = splines.eval_spline(sk, sv, knots, d=d)
        error = abs(y - values[:,d])
        print(f"f^({d}) max error: {max(error):.3e}")

    if VISUALIZE_TEST:
        print()
        print("sk: ",sk)
        print("sv: ",sv)
        print()
        padding = 0 #.1
        x = np.linspace(min(knots)-padding,max(knots)+padding,(len(knots)-1)*100+1)
        y = splines.eval_spline(sk, sv, x, d=0)

        # Make a pretty picture of the B-spline and its derivatives.
        p = Plot("Hermite Interpolant")
        p.add("Spline", x, y, mode="lines", color=p.color(1), group="s")
        p.add("Knots", knots, values[:,0], color=p.color(1), group="s")
        # Add all interesting derivatives.
        styles = ["dot", "dash", "dashdot"]
        for d in range(1, values.shape[1]):
            color = p.color(d-1) if d == 1 else p.color(d)
            dash = styles[(d-1)%len(styles)]
            dy = splines.eval_spline(sk, sv, x, d=d)
            p.add(f"D{d}", x, dy, mode="lines", color=color, dash=dash,group=d)
            p.add(f"k{d}", knots, values[:,d], color=color, group=d)
        p.show(file_name=f"spline_test-N{len(values)}-C{len(values[0])}.html")
 
    
    # 0 0 0 0 1
    # 0 0 0 1 2
    # 0 0 1 2 2
    # 0 1 2 2 2
    # 1 2 2 2 2
