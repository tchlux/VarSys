import os
from util.data import Data

all_data_file = "performance_data.pkl"
if not os.path.exists(all_data_file):
    d1 = Data.load("Pima.csv")
    d2 = Data.load("Navajo.csv")
    d1["Computer"] = ("Pima" for row in d1)
    d2["Computer"] = ("Navajo" for row in d2)
    d = d1[:]
    d += d2
    d.save("performance_data.csv")
    d.save("performance_data.pkl")

# Load data from file and reduce it to the unique tests that were run.
d = Data.load(all_data_file)
tests = d.names[:-2]
d = d[tests].unique().collect(d)
d["Num Trials"] = map(len, d["Time"])
d[d["Optimization"] == None, "Optimization"] = ""

# Print out summary information about the data.
d.summarize()
print(d)

d = d[d["Multiplicity"] <= 2]
d = d[d["Version"] == "bs-manual.f90"]
for c in sorted(set(d["Compiler"])):
    print()
    print(c)
    print(d[d["Compiler"] == c, d.names[:-4]].unique()[::-1])

exit()

# Check for nonuniform numbers of executions for different settings.
if len(set(d["Num Trials"])) > 1:
    print()
    print("RUNS ARE STILL NEEDED..")
    to_run = d[d["Num Trials"] < 20]
    num_trials = list(to_run["Num Trials"])
    to_run = to_run[to_run.names[:-3]]
    to_run["Runs Needed"] = (20-v for v in num_trials)
    # print(to_run)
    c = to_run.names.index("Compiler")
    o = to_run.names.index("Optimization")
    v = to_run.names.index("Version")
    to_run = to_run["Compiler", "Optimization", "Version"].unique().collect(to_run)
    to_run["Runs Needed"] = (sorted(set(runs))[-1] for runs in to_run["Runs Needed"])
    to_run.pop("Multiplicity")
    to_run.pop("Element")
    to_run.pop("Num Points")
    print(to_run)
    for (c,o,v,r) in to_run:
        string = "%s -lblas -llapack %s -o test real_precision.f90 test_boxspline.f90 %s"
        for i in range(r):
            print(string%(c,o,v))
    exit()

# 
from util.plot import Plot
