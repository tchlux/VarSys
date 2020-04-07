import sys, os
import numpy as np
from monotone import monotone_quintic_spline

# Set
stdout = sys.stdout
devnull = open(os.devnull,"w")

# Generate random data to construct CDFs.
num_points = 16
num_trials = 100
# 0 - 10
# 1 - 20
# 2 - 2
# 0 - 100
# Generate the random data.
np.random.seed(0)
x = np.random.random(size=(num_trials,num_points))
y = np.random.random(size=(num_trials,num_points))
# Sort the x and y data to make the sequences monotone.
x.sort(axis=1)
y.sort(axis=1)

# Set some configurataions.
max_steps = 1000
show_plot = False

def test_adaptive(adaptive=1, display=True):
    try:    adaptive = float(adaptive[0])
    except: adaptive = float(adaptive)
    # Create storage for results.
    from util.data import Data
    output_file = f"results-[adaptive-step-{adaptive}]-[{num_points}-points].csv"
    results = Data(names=["Policy", "Max Steps", "Adaptive", "Min Changes", "Max Changes", "Avg Changes", "L2"],
                   types=[str,      int,         float,      int,           int,           float,         float])
    # Load the file if it exists already (append to it).
    if os.path.exists(output_file): results = Data.load(output_file)
    if display: print()
    # Build models and measure the number of changes required to fix the data.
    for i in range(num_trials):
        if display: print(f"  i: {i:2d}", end="  --  ")
        sys.stdout = devnull
        spline, checks, changes, initial = monotone_quintic_spline(
            x[i].copy(), y[i].copy(), verbose=True,
            max_steps=max_steps, adaptive=adaptive)
        sys.stdout = stdout

        if show_plot:
            from util.plot import Plot
            p = Plot()
            p.add("points", x[i], y[i])
            p.add_func("initial", initial, [x[i][0], x[i][-1]])
            p.add_func("final  ", spline, [x[i][0], x[i][-1]])
            p.show()


        # Compute the L2 between the original spline and the final spline.
        L2 = (spline - initial)**2
        L2_error = L2.derivative(-1)(x[i][-1])
        num_changes = changes.values()
        minimum = min(num_changes)
        maximum = max(num_changes)
        average = sum(num_changes) / len(num_changes)
        results.append(["Adaptive", max_steps, adaptive, minimum, maximum, average, L2_error])
        if display: print(f" L2: {float(L2_error):.4e}  --  (min, max, avg): ({minimum:4d}, {maximum:4d}, {average:6.1f})")
    if display: print("Output in file: ",output_file)
    results.save(output_file)
    l2 = sum(results["L2"]) / len(results)
    steps = sum(results["Avg Changes"]) / len(results)
    # Return the "quality" of this "adaptive" value.
    return l2 * steps


values = [
    000.01,
    000.05,
    000.10,
    000.50,
    001.00,
    005.00,
    010.00,
    050.00,
    100.00,
    500.00
]

from util.parallel import map

run_test = lambda a: print(f"test_adaptive({a}): ",test_adaptive(a))
for _ in map(run_test, values): pass

devnull.close()
