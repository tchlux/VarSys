import numpy as np
from scipy.interpolate import PchipInterpolator

from search_monotone import monotone_quintic_spline
from fraction import Fraction
from test_points import TESTS

from util.plot import Plot
from util.system import Timer
from util.data import Data

# Import the splines Fortran code.
import og_fmodpy as fmodpy
splines = fmodpy.fimport("splines.f90", verbose=True,
                         autocompile_extra_files=True,
                         module_link_args=["-lblas", "-llapack"])


# Get the distribution of times required for doing PMQSI fits on all
# test problems, and random problems of a few different sizes.
# 
# Define an infinitely continuous "truth" for all test problems.
# 
# Test all problems with randomly spaced points and uniformly spaced
# point sets of varying sizes:
#     3  4   5   6   7    8    9    10  
#    (8, 16, 32, 64, 128, 256, 512, 1024)
# Random point sets of varying sizes:
#     3  4   5   6   7    8    9    10  
#    (8, 16, 32, 64, 128, 256, 512, 1024)
# 

trials = 100
d = Data(names=["Test Name", "N", "Trial", "Equispaced", "PMQSI fit", "PCHIP fit", "PMQSI eval", "PCHIP eval"],
         types=[str,         int, int,     bool,         float,       float,       float,        float])
t = Timer()


# p = Plot()
# test_name = "watson_test"
# N = 18
# x = np.linspace(0, 1, N)
# # Evaluate the "y" values at all points.
# y = np.asfortranarray( TESTS[test_name](x), dtype=float )
# p.add("Truth", x, y, mode="markers+lines")
# # Fit the PMQSI.
# t.start() ; sk, sc, info = splines.pmqsi(x, y) ; t.stop()
# pmqsi_fit = t.total
# # Eval the PMQSI.
# eval_pts = np.asfortranarray( np.linspace(0, 1, trials), dtype=float )
# t.start() ; splines.eval_spline(sk, sc, eval_pts, d=0) ; t.stop()
# pmqsi_eval = t.total / trials
# # Fit the PCHIP.
# t.start() ; f = PchipInterpolator(x, y) ; t.stop()
# pchip_fit = t.total
# # Eval the PCHIP.
# eval_pts = np.asfortranarray( np.linspace(0, 1, trials), dtype=float )
# t.start() ; f(eval_pts) ; t.stop()
# pchip_eval = t.total / trials


# p.add("PMQSI", np.linspace(0, 1, 1000),
#       splines.eval_spline(sk, sc, np.linspace(0, 1, 1000), d=0)[0], 
#       mode="lines")
# p.add_func("PCHIP", f, [0, 1])
# p.show()
# exit()


x = np.array([0.        , 0.00787402, 0.01574803, 0.02362205, 0.03149606,
       0.03937008, 0.04724409, 0.05511811, 0.06299213, 0.07086614,
       0.07874016, 0.08661417, 0.09448819, 0.1023622 , 0.11023622,
       0.11811024, 0.12598425, 0.13385827, 0.14173228, 0.1496063 ,
       0.15748031, 0.16535433, 0.17322835, 0.18110236, 0.18897638,
       0.19685039, 0.20472441, 0.21259843, 0.22047244, 0.22834646,
       0.23622047, 0.24409449, 0.2519685 , 0.25984252, 0.26771654,
       0.27559055, 0.28346457, 0.29133858, 0.2992126 , 0.30708661,
       0.31496063, 0.32283465, 0.33070866, 0.33858268, 0.34645669,
       0.35433071, 0.36220472, 0.37007874, 0.37795276, 0.38582677,
       0.39370079, 0.4015748 , 0.40944882, 0.41732283, 0.42519685,
       0.43307087, 0.44094488, 0.4488189 , 0.45669291, 0.46456693,
       0.47244094, 0.48031496, 0.48818898, 0.49606299, 0.50393701,
       0.51181102, 0.51968504, 0.52755906, 0.53543307, 0.54330709,
       0.5511811 , 0.55905512, 0.56692913, 0.57480315, 0.58267717,
       0.59055118, 0.5984252 , 0.60629921, 0.61417323, 0.62204724,
       0.62992126, 0.63779528, 0.64566929, 0.65354331, 0.66141732,
       0.66929134, 0.67716535, 0.68503937, 0.69291339, 0.7007874 ,
       0.70866142, 0.71653543, 0.72440945, 0.73228346, 0.74015748,
       0.7480315 , 0.75590551, 0.76377953, 0.77165354, 0.77952756,
       0.78740157, 0.79527559, 0.80314961, 0.81102362, 0.81889764,
       0.82677165, 0.83464567, 0.84251969, 0.8503937 , 0.85826772,
       0.86614173, 0.87401575, 0.88188976, 0.88976378, 0.8976378 ,
       0.90551181, 0.91338583, 0.92125984, 0.92913386, 0.93700787,
       0.94488189, 0.95275591, 0.96062992, 0.96850394, 0.97637795,
       0.98425197, 0.99212598, 1.        ])
y = np.array([0.84428838, 0.5316782 , 0.84361549, 0.3297786 , 0.35008627,
       0.70847721, 0.69155634, 0.1111936 , 0.30897958, 0.12869928,
       0.23326482, 0.37207053, 0.98028002, 0.37121875, 0.44711997,
       0.07639929, 0.37124655, 0.30264977, 0.36597435, 0.11852944,
       0.65689953, 0.07767499, 0.69625678, 0.23593911, 0.69682812,
       0.44850548, 0.96043262, 0.38980879, 0.81791332, 0.9391724 ,
       0.03734156, 0.59964689, 0.25986194, 0.2342354 , 0.51148883,
       0.94079873, 0.06192467, 0.95913663, 0.80659294, 0.90020791,
       0.05035693, 0.30073422, 0.12222783, 0.40154722, 0.17355181,
       0.06913691, 0.64696704, 0.53399559, 0.88195403, 0.1560406 ,
       0.34787674, 0.58272471, 0.08096962, 0.89225311, 0.96592443,
       0.53818353, 0.82791104, 0.49168403, 0.09566665, 0.40192077,
       0.52555024, 0.73233218, 0.13527872, 0.44023066, 0.6956519 ,
       0.63288769, 0.08403687, 0.33643686, 0.04784907, 0.14961727,
       0.20902615, 0.41518335, 0.50084075, 0.82288053, 0.77770131,
       0.49996437, 0.30945583, 0.13288597, 0.88337348, 0.45628047,
       0.61846899, 0.69548318, 0.74493713, 0.82341047, 0.49056902,
       0.43324822, 0.75290922, 0.36134999, 0.05458695, 0.81594529,
       0.65318025, 0.45910559, 0.98029373, 0.67619706, 0.13392473,
       0.35077372, 0.53499459, 0.79238128, 0.85235291, 0.59130891,
       0.0772845 , 0.37495994, 0.33805598, 0.88018508, 0.15353797,
       0.2635279 , 0.04728418, 0.46589816, 0.16144195, 0.97652471,
       0.43143275, 0.90590855, 0.20260423, 0.57158473, 0.74108906,
       0.79479433, 0.66552517, 0.99901986, 0.7162922 , 0.61176071,
       0.87479339, 0.71821783, 0.13031115, 0.76900987, 0.1506704 ,
       0.14867378, 0.29046   , 0.74657833])

print("Fitting PMQSI..")
sk, sc, info = splines.pmqsi(x, y)
exit()


# Run all tests.
for test_name in sorted(TESTS):
    print("test_name: ",test_name, " "*20)
    for N_log2 in range(3, 9):
        N = 2**N_log2
        print("  N: ",N, " "*20)
        np.random.seed(0)
        # Seed before each batch of trials.
        for trial in range(trials):
            row = [test_name, N, trial]
            for equispaced in (True, False):
                # Generate the points.
                if equispaced: x = np.linspace(0, 1, N)
                else:          x = np.random.random(size=(N,))
                x.sort()
                # Evaluate the "y" values at all points.
                y = np.asfortranarray( TESTS[test_name](x), dtype=float )
                if (test_name == "random_monotone") and (trial == 33) and \
                   (N == 128) and equispaced:
                    print()
                    print(repr(x))
                    print(repr(y))
                    print()
                # Fit the PMQSI.
                print(f"   {trial} {str(equispaced)[:1]} -- PMQSI fit", flush=True, end="\r")
                t.start() ; sk, sc, info = splines.pmqsi(x, y) ; t.stop()
                pmqsi_fit = t.total
                # Eval the PMQSI.
                print(f"   {trial} {str(equispaced)[:1]} -- PMQSI eval", flush=True, end="\r")
                eval_pts = np.asfortranarray( np.linspace(0, 1, trials), dtype=float )
                t.start() ; splines.eval_spline(sk, sc, eval_pts, d=0) ; t.stop()
                pmqsi_eval = t.total / trials
                # Fit the PCHIP.
                print(f"   {trial} {str(equispaced)[:1]} -- PCHIP fit", flush=True, end="\r")
                t.start() ; f = PchipInterpolator(x, y) ; t.stop()
                pchip_fit = t.total
                # Eval the PCHIP.
                print(f"   {trial} {str(equispaced)[:1]} -- PCHIP eval", flush=True, end="\r")
                eval_pts = np.asfortranarray( np.linspace(0, 1, trials), dtype=float )
                t.start() ; f(eval_pts) ; t.stop()
                pchip_eval = t.total / trials
                # Store these results.
                d.append( row + [equispaced, pmqsi_fit, pchip_fit, pmqsi_eval, pchip_eval] )


print(d)
d.save('temp-data.csv')
exit()

# # Get the knots and values for the test.
# knots, values = TESTS[sorted(TESTS)[-1]]
# knots = np.asfortranarray(knots, dtype=float)
# values = np.asfortranarray(values, dtype=float)
# num_knots = len(knots)

print()
print(f"knots ({len(knots)}):\n  {str(knots[:3])[:-1]} ... {str(knots[-3:])[1:]}")
print()
print(f"values ({len(values)}):\n  {str(values[:3])[:-1]} ... {str(values[-3:])[1:]}")
print()

# Compute the spline knots and spline coefficients.
t.start()
sk, sc, info = splines.pmqsi(knots, values)
t.stop()
print("info: ",info)
print(f"Fortran time: {t()}s")
print()

# Construct the PCHIP approximation.
t.start()
f = PchipInterpolator(knots, values)
t.stop()
print(f"PCHIP time: {t()}s")
print()



# # --------------------------------------------------------
# from search_monotone import monotone_quintic_spline
# t.start()
# truth = monotone_quintic_spline(knots, values, exact=True)
# t.stop() ; print(f"Python time: {t()}s")
# # --------------------------------------------------------

# g = lambda x: splines.eval_spline(sk, sc, np.array([x], order='f', dtype=float), d=0)[0]
# p = Plot()
# p.add_func("PCHIP", f, [min(knots), knots[:10][-1]])
# p.add_func("PMQSI", g, [min(knots), knots[:10][-1]])
# p.show()

# # Make a pretty picture of the B-spline and its derivatives.
# padding = 0 #.1
# x = [knots[0]]
# for i in range(len(knots) -1):
#     x += list(np.linspace(knots[i], knots[i+1], 10))[1:]
# x = np.array(x, dtype=float)
# y, info = splines.eval_spline(sk, sc, x.copy(), d=0)
# print("y: ",y)
# p = Plot("Polynomial Interpolant")
# p.add("Spline", x, y, mode="lines", color=p.color(1))
# true_y = list(map(truth, x))
# p.add("Truth", x, true_y, mode="lines", color=p.color(1), dash="dot",fill="toprevy")
# p.add("Knots", knots, values, color=p.color(1))
# p.show(file_name=f"spline_test-N{len(values)}-C{len(values)}.html")

