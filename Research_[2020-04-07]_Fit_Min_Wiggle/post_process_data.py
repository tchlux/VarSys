import os
from util.data import Data


out_file = "all_data.csv"

if os.path.exists(out_file): d = Data.load(out_file, sample=None)
else:

    # Load in all of the existing data by seed file.
    d = Data()
    for f in sorted(os.listdir()):
        if not (("seed-" == f[:len("seed-")]) and (f[-4:] == ".csv")): continue
        f_data = Data.load(f)
        f_data["Seed"] = int(f[5:-4])
        d += f_data
    # Copy in the names.
    d.reorder(d.names[:2] + ["Seed"])
    d.names[1] = "Trial Number"

    # Re-organize the accumulated seed data into one big data.
    methods = sorted({" ".join(n.split()[:2]) for n in d.names[3:]})
    names = d.names[:3] + ["Method"] + sorted(
        {" ".join(n.split()[2:]) for n in d.names[3:]} )

    # Get the indices of columns for each method.
    method_columns = {
        m : [i for i in range(len(d.names)) if m in d.names[i]]
        for m in methods }

    # Sort the indices of the method columns by the order in "names".
    for m in method_columns:
        method_columns[m].sort(key=lambda i: names.index(" ".join(d.names[i].split()[2:])))

    # Initialize the new data.
    old_d = d
    new_d = Data(names=names)
    # print("methods: ",methods)
    # print("names: ",names)
    # print("method_columns: ",method_columns)

    # Parse the old rows into new rows.
    for row in old_d:
        n, t, s = row[:3]
        for m in methods:
            new_d.append( [n,t,s,m] + [row[i] for i in method_columns[m]] )


    # Rename to the variable "d".
    d = new_d
    d.save("all_data.csv")


print(d)
methods = sorted(set(d["Method"]))
measures = d.names[4:]
max_len = max(map(len,methods))

from util.stats import rank_probability
for meth in methods:
    meth_data = list(d[d["Method"] == meth]["L2"])
    other_data = [list(d[d["Method"] == other_meth]["L2"]) for
                  other_meth in methods if other_meth != meth]
    rp = 100*rank_probability(meth_data, other_data, rank=0)
    print(f"{meth:{max_len}s} rank 0 probability: {rp:.1f}%")
exit()

import numpy as np
from util.plot import Plot
from util.stats import cdf_fit
for meas in measures:
    print(meas)
    p = Plot("",meas)
    for meth in methods:
        sub_data = d[d["Method"] == meth]
        print("",meth)
        # print(sub_data)
        print()
        fit = cdf_fit(list(sub_data[meas]))
        y = np.linspace(0,1,1000)
        x = np.percentile(list(sub_data[meas]), 100*y)
        p.add(meth, x, y, mode="lines")
    p.show(append=True)


