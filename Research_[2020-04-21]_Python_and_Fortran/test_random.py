import numpy as np
from util.plot import Plot
from search_monotone import monotone_quintic_spline


np.random.seed(2)
Ns = (4,8,16,32,64)
trials=10


N = 16
x = np.random.random(size=(N,))
y = np.random.random(size=(N,))
x.sort(); y.sort()
local_fits = {}
fit = monotone_quintic_spline(x, y, local_fits=local_fits)

print("plotting..")
p = Plot()
p.add("Points", x, y)
p.add_func("Fit", fit, [min(x), max(x)], plot_points=5000)

for i in range(len(x)):
    if i not in local_fits:
        print(f"Local fits missing {i}, skipping..")
        continue
    f, bounds = local_fits[i]
    if (i > 0): bounds[0] = (x[i-1] + x[i]) / 2
    if (i+1 < len(x)):  bounds[1] = (x[i+1] + x[i]) / 2
    p.add_func(f" local approx about {i}", f, bounds, dash="dot", opacity=0.4)

p.show()
exit()

for N in Ns:
    for t in range(trials):
        print("N, t: ",N, t)
        x = np.random.random(size=(N,))
        y = np.random.random(size=(N,))
        x.sort(); y.sort()
        failed = []
        fit = monotone_quintic_spline(x, y, failed=failed)
        if (len(failed) > 0):
            p = Plot()
            p.add("Points", x, y)
            p.add_func("Fit", fit, [min(x), max(x)])
            p.show(append=True, show=False)
# Repeate the last plot to show all of them.
p.show(append=True, show=True)


# use an arithmatic series as the step sizes.
# use SAVE attribute on things that are recalculated
# use OPTIONAL "first_call" argument to do preparations


