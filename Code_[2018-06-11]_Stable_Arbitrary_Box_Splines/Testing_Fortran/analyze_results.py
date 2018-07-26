from util.data import Data

d = Data.load("performance_data.pkl")
tests = d.names[:-2]
d = d[tests].unique().collect(d)
d["Num Trials"] = map(len, d["Time"])
d.pop("Computer")
d[d["Optimization"] == None, "Optimization"] = ""

# d.summarize()
# print(d)

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

