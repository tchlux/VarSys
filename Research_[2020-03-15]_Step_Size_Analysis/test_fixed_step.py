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

num_points = 16
num_trials = 20
# Generate the random data.
x = np.random.random(size=(num_trials,num_points))
y = np.random.random(size=(num_trials,num_points))
# Sort the x and y data to make the sequences monotone.
x.sort(axis=1)
y.sort(axis=1)

from util.data import Data
results = Data(names=["Policy", "Max Steps", "Min Changes", "Max Changes", "Avg Changes", "L2"],
               types=[str,      int,         int,           int,           float,         float])


possible_max_steps = [2**i for i in range(2,10+1)]

for max_steps in possible_max_steps:
    print("max steps:", max_steps)
    # Build models and measure the number of changes required to fix the data.
    for i in range(num_trials):
        print(f"  i: {i:2d}", end="  --  ")
        sys.stdout = devnull
        spline, checks, changes, initial = monotone_quintic_spline(
            x[i].copy(), y[i].copy(), verbose=True, max_steps=max_steps)
        sys.stdout = stdout
        # Compute the L2 between the original spline and the final spline.
        L2 = (spline - initial)**2
        L2_error = L2.derivative(-1)(x[i][-1])
        num_changes = changes.values()
        minimum = min(num_changes)
        maximum = max(num_changes)
        average = sum(num_changes) / len(num_changes)
        results.append(["Fixed Step", max_steps, minimum, maximum, average, L2_error])
        print(f" L2: {float(L2_error):.4e}  --  (min, max, avg): ({minimum:4d}, {maximum:4d}, {average:6.1f})")

print(results)
output_file = f"results-[fixed-step]-[{num_points}-points].csv"
print("Output in file: ",output_file)
results.save(output_file)

devnull.close()
