import numpy as np

# Include the "code" directory into the path and import code custom to this project.
import sys, os
sys.path += [os.path.abspath("code")]
from monotone import monotone_quintic_spline, monotone_cubic_spline
import monotone
from fraction import Fraction
from data import Data
from plot import Plot
sys.path.pop(-1)


SEED = 0
POINTS = 10

results_file = "e3_results.pkl"
results = Data(names=["num points", "seed", "total checks", "total fixes", "search steps", "checks dict", "fixes dict"],
               types=[int,          int,    int,            int,           int,            dict,          dict])

# Try to load existing results if they exist.
try: results.load(results_file)
except: pass

# Run the experiment and collect the data.
for points in range(2, 101):
    print("points: ",points)
    for seed in range(100):
        print("  seed: ",seed, end="\r")
        # Skip already checked values.
        if ((points,seed) in results["num points","seed"]): continue
        # Generate random data to test the monotonic fit function.
        np.random.seed(seed)
        x = np.linspace(0, 1, points)
        y = sorted(np.random.normal(size=(points,)))
        x -= min(x); x /= max(x)
        y -= min(y); y /= max(y); y -= max(y)
        # Convert these arrays to exact arithmetic.
        x = list(map(Fraction, x))
        y = list(map(Fraction, y))
        interval = [float(min(x)), float(max(x))]
        # Compute the monotone quintic spline.
        stdout = sys.stdout
        # Block all print statements from the function.
        with open(os.devnull, "w") as devnull:
            sys.stdout = devnull
            monotone.IS_MONOTONE_CALLS = 0
            f, checks, changes = monotone_quintic_spline(x,y=y, verbose=True)
        sys.stdout = stdout
        # Store the results.
        results.append(
            [points, seed, sum(checks.values()), sum(changes.values()),
             monotone.IS_MONOTONE_CALLS, checks, changes]
        )
    # Save the intermediate results.
    results.sort()
    results.save(results_file)

print(" "*70)

print(results[results["num points"] > 10])

d = results
d.stack(d.names[1:])
print(d)

p = Plot()
p.add("Total Checks", d["num points"], list(map(np.mean, d["total checks"])))
p.add("Total Fixes", d["num points"], list(map(np.mean, d["total fixes"])))
p.add("Search Steps", d["num points"], list(map(np.mean, d["search steps"])))
p.show()
