import os
from util.stats import cdf_fit
from util.plot import Plot
from util.data import Data



CREATE_ALL_DATA_PLOT = True
if CREATE_ALL_DATA_PLOT:
    # Otherwise load the combined data file.
    d = Data.load("combined-data.csv", sample=False)
    print("Reducing data..", flush=True)
    # d = d[d["host"] == "m03"].copy()
    d = d[d["freq"] == 1600000].copy()
    d = d[d["fsize"] == 1024].copy()
    d = d[d["rsize"] == 4].copy()
    d = d[d["threads"] == 32].copy()
    d = d[d["test"] == "readers"].copy()
    print(d.summary)

    print("Making plot..")
    p = Plot("Various Hosts, Freq 1.6, Fsize 1024, Rsize 4, Threads 32, Readers",
             "Read Throughput", "CDF Value")
    for data_name in sorted(set(d["data"])):
        v = d[d['data'] == data_name]
        fit = cdf_fit(v["throughput"])
        p.add_func(data_name, fit, fit())
    p.show(file_name="[2019-04-27]_data_set_comparison.html")


BUILD_ALL_DATA = False
if BUILD_ALL_DATA:
    if (not os.path.exists("combined-data.pkl")):
        names = ["old-cpu-scaling.csv", "li-full.csv",
                 "rt-kernel.csv", "nort-nojournal.csv", "nort-perf.csv"]
        keep_cols = ["host", "media", "freq", "fsize", "rsize", "threads", "test", 
                     "iter", "throughput", "runtime"]
        config_cols = keep_cols[:7]
        unique_configs = None
        all_data = Data(names=["data"] + keep_cols)

        for n in names:
            print("\n"*3)
            print(n)
            try:
                d = Data.load(n+".pkl")
            except: 
                d = Data.load(n, sample=None)
                # Remove the index column.
                if "" in d.names: d.pop("")
                # d.stack([n for n in d.names if n not in config_cols])
                d.save(n+".pkl")

            # Reduce this particular data set to only one host.
            if (n == "li-full.csv"): d = d[d["host"] == "m03"].copy()

            # If the unique configurations have not been identified, find them.
            if (type(unique_configs) == type(None)):
                print("Collecting unique configurations..")
                unique_configs = set(map(tuple, d[config_cols[2:]]))
            else:
                print("Reducing to desired configs..")
                d = d[(i for i in range(len(d)) if
                       tuple(d[i,config_cols[2:]][0]) in unique_configs)]

            d.stack([n for n in d.names if n not in config_cols])
            d = d[keep_cols].copy()
            # Add this data to the master data.
            for row in d: all_data.append( [n[:-4]] + list(row) )

        print("Loading old data..", flush=True)
        d = Data.load("2018-stacked.pkl")
        print("Reducing old data..", flush=True)
        d = d[(i for i in range(len(d)) if
               tuple(d[i,config_cols[2:]][0]) in unique_configs)]
        for row in d: all_data.append( ["2018-data"] + list(row) )

        all_data.save("combined-data.pkl")
        d = all_data
        stacked = [n for (n,t) in zip(d.names, d.types) if t == list]
        print("Unstacking combined data..")
        d.unstack(stacked)
        d.save("combined-data.csv")
        print(d)


EXAMINE_LI_FULL = False
if EXAMINE_LI_FULL:
    n = "li-full.csv"
    try: d = Data.load(n+".pkl")
    except: 
        d = Data.load(n)
        d.save(n+".pkl")

    # Test for effects, just to look at what is different.
    d.effect("host")# , use_ks=True)
    # Plot the throughput distributions
    hosts = [(h, list(d[d["host"]==h]["throughput"]))
             for h in sorted(set(d["host"]))]
    p = Plot("throughput CDFs by host")
    for h, dist in hosts:
        f = cdf_fit(dist)
        p.add_func(h, f, f())
    p.show(file_name="[2019-04-03]_li-full_analysis.html")
    # Plot the runtime distributions
    runt = [(h, list(d[d["host"]==h]["runtime"]))
            for h in sorted(set(d["host"]))]
    p = Plot("runtime CDFs by host")
    for h, dist in runt:
        f = cdf_fit(dist)
        p.add_func(h, f, f())
    p.show(append=True)
