import os, pickle
import numpy as np
from util.data import read_struct
from util.plotly import Plot
from util.stats import cdf_fit_func

CLEAN_DATA_CSV = "clean_data.csv"
CLEAN_DATA_CSV_PKL = "clean_data.csv.pkl"

if not os.path.exists(CLEAN_DATA_CSV_PKL):
    data = read_struct(CLEAN_DATA_CSV)
    print("Saving '%s'..."%(CLEAN_DATA_CSV_PKL))
    with open(CLEAN_DATA_CSV_PKL, "wb") as f:
        pickle.dump(data,f)
else:
    print("Loading '%s'..."%(CLEAN_DATA_CSV_PKL))
    with open(CLEAN_DATA_CSV_PKL, "rb") as f: 
        data = pickle.load(f)

header = list(data.dtype.names)
h_config = header[:4]
h_test = header[4]
h_values = header[5:]
tests = sorted(np.unique(data[h_test]))
test_len = max(map(len, tests))
# Extract a subset of data, only the "readers" test
# data = data[data["Test"] == "readers"]

# ============================================================

count = 0
failures = 0
slow_configs = []
throughput_threshold = 40000
total_threshold = 80

p = Plot("Low throughput configurations", "Sample number", "Throughput")
for row in data:
    total = sum(1 for h in h_values if row[h] < throughput_threshold)
    if total > 0:
        slow_configs.append((tuple(row[h] for h in h_config), total))
        if total >= total_threshold:
            raw_values = np.array([row[h] for h in h_values])
            name = str(slow_configs[-1][0] + (row["Test"],))
            # p.add("",np.arange(len(raw_values)), raw_values, 
            #       frame=name, group="")
            p.add_histogram("", raw_values, frame=name, group="")
    if (row["File Size"] == 4 and row["Record Size"] == 4): 
        count += 1
        if total > 0:
            failures += total
slow_configs = sorted(slow_configs,key=lambda r: r[1])[-16:]
for (c,t) in slow_configs:
    print(c,t)

p.plot(show_slider_labels=False, frame_label="System Configuration: ",
       data_easing=False, show_legend=False, bounce=False,
       redraw=True, show_play_pause=False)
exit()

print()
print("Number of slow runs: ", failures)
print("Number of (4,4) runs:", count*150)
print("Total number of runs:", len(data)*150)
print()
all_throughputs = (np.vstack([data[h] for h in h_values]).T).flatten()

print("Total tests:        ",len(all_throughputs))
print("Smallest throughput:",min(all_throughputs))
print("Throughputs <=%s: "%throughput_threshold,
      sum(np.where(all_throughputs <= 10000, 1, 0)))

exit()
# ============================================================
