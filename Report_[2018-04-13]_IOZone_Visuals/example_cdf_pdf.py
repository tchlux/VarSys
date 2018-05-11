import os
import numpy as np
from util.data import Struct
from util.plotly import Plot, multiplot
from util.stats import cdf_fit_func

test = "initial_writers"

# Brightness and color series
# 
# 0 -- "Frequency"      "initial_writers"
# 1 -- "File Size"      "rewriters"      
# 2 -- "Record Size"    "readers"        
# 3 -- "Num Threads"    "re-readers"     
#                       "random_readers" 
#                       "random_writers" 


# color_by = 1
# brightness_by = 2
# reduction = 2800000
# reduction_by = "Frequency"

color_by = 3
brightness_by = 2
reduction = 1048576
reduction_by = "File Size"


######################################################################
if not os.path.exists(f"data_{test}_[cdfs].dill.gz"):
    # Otherwise, load the pre-computed data file.
    print("Loading full data...")
    data = Struct().load("data_full.pkl.gz")
    print("Reducing only to columns of interest...")
    data = data[[n for i,n in enumerate(data.names) if i in {2,3,4,5,7,8}]]
    print(f"Reducing to rows with test \"{test}\"...")
    data = data[data["Test"] == test]
    print("Collecting repeated trials...")
    counts = {}
    for row in data:
        setting = tuple(row[:-1])
        counts[setting] = counts.get(setting, []) + [row[-1]]
    print()
    print("Unique configs: ", len(counts))
    print("Unique num runs:", sorted(set(len(counts[s]) for s in counts)))
    print("Less than 150:  ", sum(1 for s in counts if len(counts[s]) < 150))
    print()
    # Remove configurations with too little data
    for s in list(counts):
        if len(counts[s]) < 150: 
            print(f"Too little data ({counts[s]}):", s)
            counts.pop(s)
        else:
            counts[s] = counts[s][:150]

    # Process raw data into CDF functions.
    print()
    print("Generating CDFs for each configuration.")
    dists = Struct(names=(data.names[:-1] + ["Throughput CDF"]))
    for s in sorted(counts):
        dists.append( list(s) + [cdf_fit_func(counts[s])] )
    print("Done.")
    print()
    dists.save(f"data_{test}_[cdfs].dill.gz")
else:
    dists = Struct().load(f"data_{test}_[cdfs].dill.gz")
######################################################################


np.random.seed(0)
ind = np.random.randint(len(dists))
f = dists[ind][-1]
p1 = Plot("","File I/O Write Throughput (kb/s)","Cumulative Distribution Function Value")
p1.add_func("", f, f(), color=p1.color(1))
p2 = Plot("","File I/O Write Throughput (kb/s)","Percentage of Runs")
step = f()[1] - f()[0]
j = .1
spots = lambda x: np.linspace(x-j*step, x+j*step, 1000)
deriv = lambda x: (f(spots(x)+.01*step).mean() - f(spots(x)-.01*step)).mean() / (.02*step)
p2.add_func("", lambda x: 100000*deriv(x), f(), fill="tozeroy")
multiplot([[p1],[p2]], show_legend=False, shared_x=True, file_name="CDF_and_PDF_[2800000,1024,128,32,initial_writers].html")
print(dists[ind])

