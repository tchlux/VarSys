import os, math
from util.plot import Plot
from util.data import Data
from util.stats import cdf_fit

results_folder = os.path.join(os.path.dirname(os.path.abspath(__file__)), "Predictions")
for f in os.listdir(results_folder):
    if ("results" != f[:len("results")]): continue
    path = os.path.join(results_folder, f)
    d = Data.load(path)
    print(f, d.shape)
    # Get the target column and the list of models that made predictions.
    target_index = d.names.index("Indices") - 1
    target = d.names[target_index]
    models = [c.split()[0] for c in d.names[target_index+2:] if target in c][::-1]
    # Generate a plot of the model errors.
    p = Plot(f"{f} Absolute Error Distributions")
    max_error = -float('inf')
    for m_name in models:
        errors = [truth - guess
                  for (truth,guess) in zip(d[target], d[f"{m_name} {target}"])
                  if ((truth != None) and (guess != None))]
        errors = list(map(abs, errors))
        max_error = max(max_error, max(errors))
        p.add_box(m_name, errors)
        # print(f"{m_name} absolute relative error (min,max):", (min(errors), max(errors)))
        # fit = cdf_fit(errors)
        # p.add_func(m_name, fit, fit())
    p.show(append=True, y_range=[math.log(.3),math.log(max_error)*.5],
           y_axis_settings=dict(type="log"))
    print()

# FOR EACH DATA SET

# Absolute error box plot
# Absolute error quartiles table with bold minimum values
# Stacked bar chart per-point fit & predict times


# Absolute error distribution [(model) x (data set)].
# Absolute error distribution versus distance to nearest neighbor.
# Absolute error distribution versus diameter of influencers.
# Histogram of number of contributors for each and model.
# Table of fit times for each model.

# Box plot distribution of prediction time by dimension (and by number points)
