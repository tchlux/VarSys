from util.data import Data

file_name = "[2018-11]_tensorflow-runtime-luna.csv"

d = Data.load(file_name, sep=",", verbose=True)

# Get the columns for data and times.
normal_cols = [i for (i,n) in enumerate(d.names) if "time" not in n.lower()]
time_cols = [i for (i,n) in enumerate(d.names) if "time" in n.lower()]
# Create a single time column.
d["Time"] = ([row[i] for i in time_cols if row[i] != None] for row in d)
# Reduce data to the configs and time list.
d = d[:20, normal_cols+[-1]]
d["Samples"] = (len(t) for t in d["Time"])
d["Mean Time"] = (sum(t)/len(t) for t in d["Time"])
d["Min Time"] = (min(t) for t in d["Time"])
d["Max Time"] = (max(t) for t in d["Time"])
d.reorder(d.names[:-5] + d.names[-4:] + [d.names[-5]])
d._max_display = 100
d.sort()
print(d)
d.effect("Mean Time")
print()




from util.plot import Plot
from util.stats import cdf_fit

group_by = [None] + d.names[:-4]
for g in group_by:
    if g != None: name = "Grouped by "+str(g)
    else:         name = ""
    p = Plot(name, "Runtime", "CDF(x)")
    groups = []
    for vals in d:
        net_name, steps, batch_size, data_name, db, samples, mean_time, min_time, max_time, times = vals
        if g != None:
            group = vals[d.names.index(g)]
            config_name = f"{g}: '{group}'"
            if group not in groups: groups += [group]
            c = p.color(groups.index(group))
        else:
            config_name = f"{net_name}-{data_name}-{db} ({batch_size},{samples})"
            group = None
            c = None
        cdf = cdf_fit(times)
        p.add_func(config_name, cdf, cdf(), group=group, color=c)
    p.show(append=(g != group_by[0]),
           show=(g == group_by[-1]),
           file_name="exploring_groups.html")


# from pandas import DataFrame
# d2 = DataFrame(list(d), columns=d.names)
# print(d2.corr())

#   
# d.diff()
