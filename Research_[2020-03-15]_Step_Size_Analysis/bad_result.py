# Set up an experiment that can vary the "step size" function parameters.
#  Returns:
#    - distance from start point * number of steps needed to obtain solution


# Generate random data to construct CDFs.

import sys, os
import numpy as np
from monotone import monotone_quintic_spline

stdout = sys.stdout
devnull = open(os.devnull,"w")
np.random.seed(0)

num_points = 10
num_trials = 100
# Generate the random data.
x = np.random.random(size=(num_trials,num_points))
y = np.random.random(size=(num_trials,num_points))
# Sort the x and y data to make the sequences monotone.
x.sort(axis=1)
y.sort(axis=1)

# Build models and measure the number of changes required to fix the data.
for i in range(num_trials):
    print()
    print("i: ",i)
    # sys.stdout = devnull
    spline, checks, changes, initial = monotone_quintic_spline(
        x[i], y[i], verbose=True, max_steps=100)
    # sys.stdout = stdout
    # Convert the change in values into a numpy array.
    L2 = (spline - initial)**2

    # # Convert all functions to inexact arithmetic.
    # for f in L2.functions:
    #     f.coefficients = [float(c) for c in f.coefficients]
    # print("L2: ",L2)
    from util.plot import Plot
    p = Plot()
    min_max = x[0][0], x[0][-1]
    p.add("Points", x[0], y[0])
    p.add_func("Original", initial, min_max)
    p.add_func("Final", spline, min_max)
    p.add_func("Diff", (spline - initial), min_max)
    p.add_func("L2", L2, min_max)
    p.add_func("L2 integral", L2.derivative(-1), min_max)
    p.show()
    # 10   -> 0.08964440432913916
    # 100  -> 0.10453634164051213
    # 1000 -> 0.10453634279381023
    print("L2.derivative(-1)(x[i][-1]): ",float(L2.derivative(-1)(x[i][-1])))
    exit()

devnull.close()
