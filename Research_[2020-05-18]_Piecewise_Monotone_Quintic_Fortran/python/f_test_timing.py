import numpy as np
from scipy.interpolate import PchipInterpolator

from mqsi import monotone_quintic_spline
from test_points import TESTS

from util.system import Timer
from util.data import Data
from util.plot import Plot

SHOW_PCHIP = False
VERSION = 11.1
TRIALS = 100

# Import the splines Fortran code.
import fmodpy
spline = fmodpy.fimport("SPLINE.f90")
mqsi = fmodpy.fimport("MQSI.f90")
# Define an object that has all of the expected functions.
class splines:
    mqsi = mqsi.mqsi
    eval_spline = spline.eval_spline

# Get the distribution of times required for doing MQSI fits on all
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
# 
d = Data(names=["Test Name", "N", "Trial", "Equispaced", "MQSI fit", "PCHIP fit", "MQSI eval", "PCHIP eval"],
         types=[str,         int, int,     bool,         float,       float,       float,        float])
t = Timer()

file_name = f"v{VERSION}.csv"
import os

if not os.path.exists(file_name):
    # Run all tests.
    for test_name in sorted(TESTS):
        print("test_name: ",test_name, " "*20)
        for N_log2 in range(3, 9):
            N = 2**N_log2
            print("  N: ",N, " "*20)
            np.random.seed(0)
            # Seed before each batch of trials.
            for trial in range(TRIALS):
                row = [test_name, N, trial]
                for equispaced in (True, False):
                    # Generate the points.
                    if equispaced: x = np.linspace(0, 1, N)
                    else:          x = np.random.random(size=(N,))
                    x.sort()
                    # Evaluate the "y" values at all points.
                    y = np.asfortranarray( TESTS[test_name](x), dtype=float )
                    # Fit the MQSI.
                    print(f"   {trial} {str(equispaced)[:1]} -- MQSI fit", flush=True, end="\r")
                    nd = len(x)
                    sk = np.ones(3*nd + 6)
                    sc = np.ones(3*nd)
                    t.start() ; sk, sc, info = splines.mqsi(x, y, sk, sc) ; t.stop()
                    if (info != 0):
                        print()
                        print()
                        print("test_name:   ",test_name)
                        print("N:           ",N)
                        print("trial:       ",trial)
                        print("equispaced:  ",equispaced)
                        print("x =", repr(x))
                        print("y =", repr(y))
                        print()
                        raise(Exception(f"PQMSI FIT ERROR: {info}"))
                    mqsi_fit = t.total
                    # Eval the MQSI.
                    print(f"   {trial} {str(equispaced)[:1]} -- MQSI eval", flush=True, end="\r")
                    eval_pts = np.asfortranarray(np.linspace(min(x), max(x), TRIALS), dtype=float)
                    t.start() ; _, info = splines.eval_spline(sk, sc, eval_pts, d=0) ; t.stop()
                    if (info != 0):
                        print()
                        print()
                        print("test_name:   ",test_name)
                        print("N:           ",N)
                        print("trial:       ",trial)
                        print("equispaced:  ",equispaced)
                        print("x =", repr(x))
                        print("y =", repr(y))
                        print()
                        raise(Exception(f"PQMSI EVAL ERROR: {info}"))
                    mqsi_eval = t.total / TRIALS
                    # Fit the PCHIP.
                    print(f"   {trial} {str(equispaced)[:1]} -- PCHIP fit", flush=True, end="\r")
                    t.start() ; f = PchipInterpolator(x, y) ; t.stop()
                    pchip_fit = t.total
                    # Eval the PCHIP.
                    print(f"   {trial} {str(equispaced)[:1]} -- PCHIP eval", flush=True, end="\r")
                    eval_pts = np.asfortranarray(np.linspace(min(x), max(x), TRIALS), dtype=float)
                    t.start() ; f(eval_pts) ; t.stop()
                    pchip_eval = t.total / TRIALS
                    # Store these results.
                    d.append( row + [equispaced, mqsi_fit, pchip_fit, mqsi_eval, pchip_eval] )
    print(" "*40)
    d.save(file_name)
else:
    d.load(file_name)

print(d)

# Plot the distribution of times required to construct fits and make evaluations.
from util.stats import cdf_fit
p = Plot("Time (seconds)", "CDF")

# Get the list of versions that are available.
versions = [f for f in os.listdir() if ("v" == f[0]) and (".csv" == f[-4:])]
versions = sorted(versions, key=lambda v: float(v[1:-4]))

# --------------------------------------------------------------------
# Add MQSI approximation.
cdf = cdf_fit(d["MQSI fit"])
p.add_func(f"v{VERSION} MQSI fit", cdf, cdf(), color=1, group='current')
cdf = cdf_fit(d["MQSI eval"])
p.add_func(f"  v{VERSION} MQSI eval", cdf, cdf(), color=1, group='current')

for step, v in reversed(list(enumerate(versions))):
    # Load old data.
    d2 = Data.load(v)
    # Only add *old* MQSI information.
    if (float(v[1:-4]) == VERSION): continue
    if (int(float(v[1:-4])) == int(VERSION)): d = None
    else:                                     d = "dot"
    # Add old MQSI information.
    if ("MQSI fit" in d2): cdf = cdf_fit(d2["MQSI fit"])
    else:                  cdf = cdf_fit(d2["PMQSI fit"])
    p.add_func(f"V{v[1:-4]} MQSI fit", cdf, cdf(), color=1,
               opacity=(step+1) / (len(versions)+1), dash=d,
               group=v[1:-4].split('.')[0]+" fit")
    if ("MQSI eval" in d2): cdf = cdf_fit(d2["MQSI eval"])
    else:                   cdf = cdf_fit(d2["PMQSI eval"])
    p.add_func(f"V{v[1:-4]} MQSI eval", cdf, cdf(), color=1,
               opacity=(step+1) / (len(versions)+1), dash=d,
               group=v[1:-4].split('.')[0]+" eval")


if SHOW_PCHIP:
    # Add all the PCHIP information for comparison.
    for step, v in reversed(list(enumerate(versions))):
        # Load old data.
        d2 = Data.load(v)
        # Add old PCHIP information.
        cdf = cdf_fit(d2["PCHIP fit"])
        p.add_func(f"PCHIP fit", cdf, cdf(), color=0,
                   opacity=(step+1) / (len(versions)+1),
                   group='PCHIP fit')
        cdf = cdf_fit(d2["PCHIP eval"])
        p.add_func(f"PCHIP eval", cdf, cdf(), color=0,
                   opacity=(step+1) / (len(versions)+1),
                   group='PCHIP eval')


p.show()