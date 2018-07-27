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
    d.save(all_data_file)

# Load data from file and reduce it to the unique tests that were run.
d = Data.load(all_data_file)
tests = d.names[:-2]
d = d[tests].unique().collect(d)
d["Num Trials"] = map(len, d["Time"])
d[d["Optimization"] == None, "Optimization"] = ""
d["Computer"] = (c[0] for c in d["Computer"])
d["Time"] = d.pop("Time")

# Print out summary information about the data.
d.summarize()
print(d)
print()

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
from util.stats import cdf_fit_func
from util.plot import Plot, multiplot
import numpy as np

legend_settings = dict(
    xanchor = "center",
    yanchor = "top",
    x = .5,
    y = -.18,
    orientation = "h",
)

# Break up the data by compiler
ifort = d[d["Compiler"] == "ifort"]
gfort = d[d["Compiler"] == "gfortran"]
sun95 = d[d["Compiler"] == "/home/f/ltw/bin/ftn95.sun"]
# Set the compiler names
ifort.compiler = "ifort"
gfort.compiler = "gfortran"
sun95.compiler = "sun95"

# Order the different comparitors by which is best.
compilers = [gfort, ifort, sun95]
optimizations = sorted(set(d["Optimization"]))
versions = ["bs-dynamic.f90", "bs-automatic.f90", "bs-manual.f90"]


# Plot by optimization level for each compiler.
plots = []
for c in compilers:
    p = Plot(f"Aggregate Runtime CDFs for Different Optimizations on Each Compiler for {len(c)//4} Tests", 
             "Execution Time (sec)", c.compiler)
    for i,opt in enumerate(optimizations):
        perf_data = [v for t in c[c["Optimization"]==opt, "Time"] for v in t[0]]
        cdf = cdf_fit_func(perf_data)
        if (opt.strip() == ""): opt = "No Optimization"
        p.add_function(opt, cdf, cdf(), group=opt, show_in_legend=(c==compilers[0]))
    plots.append(p.plot(html=False, y_range=[.5,1], 
                        legend=legend_settings))
multiplot([[p] for p in plots], show=False, append=True, shared_x=True)

# 
plots = []
for c in compilers:
    name = c.compiler
    # Plot by code version for Sun.
    c = c[ c["Optimization"] == "-O3" ]
    c = c[ c["Num Points"] == "4K" ]
    first = (name==compilers[0].compiler)
    p = Plot(f"'O3' Runtime CDFs for Different Box-Spline Versions on Each Compiler<br>over {len(c)//3} Tests, Each with 4K Evaluation Points", 
             "Execution Time (sec)", f"{name}")
    for i,version in enumerate(versions):
        perf_data = [v for t in c[c["Version"]==version, "Time"] for v in t[0]]
        cdf = cdf_fit_func(perf_data)
        if ("manual" in version): version = version.replace("manual", "allocate")
        version = version.replace("bs-","").replace(".f90","")
        p.add_function(version, cdf, cdf(), color=p.color(i), group=version,
                       show_in_legend=first)
    plots.append( p.plot(html=False, legend=legend_settings, y_range=[.4,1]) )
multiplot([[p] for p in plots], append=True, shared_x=True)

