import os
from util.data import Data
from util.stats import cdf_fit
from util.plot import Plot


if not os.path.exists("delta.csv"):
    data_files = ["delta-noperf.csv", "delta-noperf-nojournal.csv",
                  "delta-perf-nojournal.csv", "delta-perf.csv"]
    keep_cols = None
    for f in data_files:
        print(f)
        d = Data.load(f, sample=None)
        d['data'] = f[:-4]
        if (type(keep_cols) == type(None)):
            keep_cols = ["","data"] + d.names[1:-1]
            merged_data = d[keep_cols].copy()
        else:
            merged_data += d
        print(d)
        print("\n"*2)
    print()
    print(merged_data)
    d = merged_data
    d.save("delta.csv")
    d.save("delta.pkl")    
    # Create a stacked version as well.
    config_cols = ['data', 'host', 'media', 'freq', 'fsize', 'rsize', 'threads', 'test']
    to_stack = [n for n in d.names if n not in config_cols]
    d.stack(to_stack)
    print("Removing empty stacks..")
    for i in range(d.shape[0]):
        for j in range(d.shape[1]):
            if (type(d[i,j]) == list) and (set(d[i,j]) == {None}):
                d[i,j] = None
    d.save("delta-stacked.pkl")
else:
    print("Loading data from file..", flush=True)
    # d = Data.load("delta.csv")
    # d = Data.load("delta.pkl")
    d = Data.load("delta-stacked.pkl")


# Pop out columns I'm not going to use.
d.pop("media")
d.pop("")
d.pop("subiter")

print()
print(d)
print()

tests = sorted(set(d["test"]))
datas = sorted(set(d["data"]))
hosts = sorted(set(d["host"]))
run_config_cols = ["freq", "fsize", "rsize", "threads"]

print("datas: ",datas)
print("hosts: ",hosts)
print("tests: ",tests)


for t in tests:
    print('-'*70)
    print("Running for test", t)
    print()
    # Reduce to the specified test.
    d_test = d[d["test"] == t]

    for n in datas:
        print()
        print("n: ",n)
        view = d_test[d_test["data"] == n]
        configs = view[run_config_cols].unique()
        configs.sort()
        p = Plot(f"{t} Througphut for {n}", "Througphut", "CDF Value")
        for config in configs:
            matching_rows = [i for i in range(len(view))
                             if all(view[i,run_config_cols] == config)]
            matching = view[matching_rows]
            machines = sorted(set(matching["host"]))
            for mach in machines:
                # Add to the plot the distributions of Throughput for this machine
                row = matching[ matching["host"] == mach ][0]
                f = cdf_fit(row["throughput"])
                p.add_func(mach, f, f(), frame=str(tuple(config)))
        p.show(frame_label="(freq, fsize, rsize, threads) -- ",
               show_play_pause=False, append=(n != datas[0]),
               file_name=f"[2019-04-29]_Data_{t}_Througput_Overview.html")
