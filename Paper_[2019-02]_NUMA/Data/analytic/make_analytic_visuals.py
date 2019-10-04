from util.data import Data

d = Data.load("2-analytic-results.pkl")

# # Strip out the invalid tests (bad number of train points).
# to_keep = [i for i in range(len(d)) if ((d[i,"Train"]-1)%10)]
# # Make a copy so that we can do assignment etc.
# d = d[to_keep][:]
# d.save("1-analytic-results.pkl")

# Generate interesting extra columns.
d["Abs Errors"] = ([float(abs(v)) for v in l] for l in d["Errors"])
d["Mean Abs Error"] = (sum(l) / len(l) for l in d["Abs Errors"])
d["Min Abs Error"] = (min(l) for l in d["Abs Errors"])
d["Max Abs Error"] = (max(l) for l in d["Abs Errors"])
# d.summarize()

config_cols = ["SNR", "Train"]

configs = d[:,config_cols].unique()
configs.sort()
configs.max_display=float('inf')
print(configs)

from util.plot import Plot
from util.stats import cdf_fit

for conf in configs:
    interest = d[d[:,config_cols] == conf]
    interest.sort()
    print(interest)
    S, N = conf
    p = Plot(f"{N} train with SNR {S}")
    seen = {"NearestNeighbor"}
    for row in interest:
        name = row["Algorithm"]
        # Skip algorithms that have already been plotted.
        if name in seen: 
            print(f"  skipping extra '{name}'")
            continue
        else:            seen.add(name)
        # Plot the performance of this algorithm.
        fit_time = row["Fit Time"]
        pred_time = row["Predict Time"]
        errors = row["Abs Errors"]
        # f = cdf_fit(errors)
        # p.add_func(name, f, f())
        p.add_box(name, errors)
    p.show(append=True, show=all(conf == configs[-1]))
