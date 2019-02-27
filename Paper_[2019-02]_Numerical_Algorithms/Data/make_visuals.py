import os, math
import numpy as np
from util.plot import Plot
from util.data import Data
from util.stats import cdf_fit
from util.misc.paper import latex_table

results_folder = os.path.join(os.path.dirname(os.path.abspath(__file__)), "Predictions")
percentiles = [0, 25, 50, 75, 100]

for f in os.listdir(results_folder):
    if ("results" != f[:len("results")]): continue
    path = os.path.join(results_folder, f)
    d = Data.load(path)
    print()
    print(f, d.shape)
    print(d)
    print()
    # Get the target column and the list of models that made predictions.
    target_index = d.names.index("Indices") - 1
    target = d.names[target_index]
    models = [c.split()[0] for c in d.names[target_index+2:] if target in c][::-1]
    # Generate a plot of the model errors.
    p = Plot(f"{f} Absolute Error Distributions")
    max_error = -float('inf')
    table = [["Algorithm", "Min"] + [f"${v}^{{th}}$" for v in percentiles[1:-1]] + ["Max"]]
    for m_name in models:
        errors = [truth - guess
                  for (truth,guess) in zip(d[target], d[f"{m_name} {target}"])
                  if ((truth != None) and (guess != None))]
        print()
        # Skip bad models.
        if 'nan' in map(str,errors): errors = [v for v in errors if str(v) != 'nan']
        if (len(errors) == 0):
            print(f"Skipping {m_name} because no predictions are found..")
            continue
        abs_errors = list(map(abs,errors))
        max_error = max(max_error, max(abs_errors))
        print("errors: ",len(errors), errors[:100])
        print("m_name: ",m_name)
        print("max_error: ",max_error)
        p.add_box(m_name, abs_errors)
        table.append([m_name] + [np.percentile(abs_errors, p) for p in percentiles])
        # print(f"{m_name} absolute relative error (min,max):", (min(errors), max(errors)))
        # fit = cdf_fit(errors)
        # p.add_func(m_name, fit, fit())
    # Fix the table to have bold entries for the minimal values at each percentile.
    col_mins = [min(row[i] for row in table[1:]) for i in range(1+len(percentiles))]
    for row in range(1,len(table)):
        for col in range(1,len(table[row])):
            is_min = (table[row][col] == col_mins[col])
            table[row][col] = f"{table[row][col]:.3e}"
            if is_min: table[row][col] = "\\textbf{" + table[row][col] + "}"
    print()
    latex_table(table[1:], table[0])
    print()
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
