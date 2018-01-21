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
            p.add("",np.arange(len(raw_values)), raw_values, 
                  frame=name, group="")
    if (row["File Size"] == 4 and row["Record Size"] == 4): 
        count += 1
        if total > 0:
            failures += total
slow_configs = sorted(slow_configs,key=lambda r: r[1])[-16:]
for (c,t) in slow_configs:
    print(c,t)

p.plot(show_slider_labels=False, frame_label="System Configuration: ",
       data_easing=True, show_legend=False, bounce=True)

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


# Generate the sorted list of indices (sorted by min value, then by range)
indices = np.arange(len(data))

PLOT_RAW_THROUGHPUT = True
SORT_BY_RANGE = True
SORT_BY_CONFIG = False

if SORT_BY_RANGE:
    sorted_by = "Range"
    indices = sorted(indices, key=lambda i: 
                     min([data[i][h] for h in h_values]))
    indices = sorted(indices,key=lambda i:
                     max([data[i][h] for h in h_values]) -
                     min([data[i][h] for h in h_values]))
elif SORT_BY_CONFIG:
    sorted_by = "Config"
    order = h_config[::-1]
    print("Using order:", order)
    for h in order:
        indices = sorted(indices, key=lambda i: data[i][h])

# Convert sorted indices back to numpy array, randomly keep a subset
indices = np.array(indices)
to_keep = np.arange(len(data))
np.random.seed(0)
np.random.shuffle(to_keep)
to_keep = np.array(sorted(to_keep[:300]))
indices = indices[to_keep]
# Create the plot

if PLOT_RAW_THROUGHPUT:
    p = Plot("Raw Throughput Values for 300 System Configurations ('Readers' Test) Sorted by %s"%sorted_by+
             "<br><i style='font-size:80%;'>Use 'Autoscale' to zoom in on a series. Use 'Home' to reset to full range. (top right button panel)"+
             "<br>Move slider by clicking and dragging to manually explore. (below system configuration text)</i>",
             "Sample Number", "Throughput")
else:
    p = Plot("CDFs of 300 System Configurations ('Readers' Test) Sorted by %s"%sorted_by+
             "<br><i style='font-size:80%;'>Use 'Autoscale' to zoom in on a CDF. Use 'Home' to reset to full range. (top right button panel)"+
             "<br>Move slider by clicking and dragging to manually explore. (below system configuration text)</i>",
             "Throughput", "CDF Value (P[Trial Throughput < Throughput])")
color = p.color()
global_min_max = [float("inf"), -float("inf")]
for i in indices:
    row = data[i]
    values = [row[h] for h in h_values]
    min_max = (min(values), max(values))
    global_min_max[0] = min(min_max[0], global_min_max[0])
    global_min_max[1] = max(min_max[1], global_min_max[1])
    throughput_range = min_max[1] - min_max[0]
    config = [h+" "+str(row[h]) for h in h_config]
    name = "%.1e"%throughput_range+("&nbsp;"*4)+"<b>System:</b> "+(("&nbsp;"*4).join(config))
    if PLOT_RAW_THROUGHPUT:
        p.add("", np.arange(len(values)), values, group="", frame=name, show_in_legend=False)
    else:
        p.add_func("", cdf_fit_func(values), min_max, group="",
                   plot_points=500, frame=name,
                   show_in_legend=False)

if PLOT_RAW_THROUGHPUT:
    p.plot(show_legend=False, frame_label="<b>Throughput Range:</b> ",
           file_name="[2018-01-15]_System_Raw_Throughput_Sample_300(0)(%s).html"%sorted_by,
           show_slider_labels=False, show_play_pause=False, initial_frame=name)
else:
    p.plot(show_legend=False, frame_label="<b>Throughput Range:</b> ",
           file_name="[2018-01-15]_System_Throughput_CDFs_Sample_300(0)(%s).html"%sorted_by,
           show_slider_labels=False, show_play_pause=False, initial_frame=name)

exit()

# Reduce to the "readers" test
# Get a list of all unique system configurations (excluding num threads)
# Sort that list of configurations by the average (max - min) range over threads
# Randomly sample down to 



# data = data[data["File Size"] == ]
# threads = sorted(set(data["Num Threads"]))
# print(threads)
# min_pts = []
# max_pts = []
# var_pts = []
# for t in threads:
#     subdata = data[data["Num Threads"] == t]
#     all_throughputs = (np.vstack([subdata[h] for h in h_values]).T).flatten()
#     print(t,min(all_throughputs), max(all_throughputs),sep="\t")
#     var_pts.append((t,np.var(all_throughputs)))
#     min_pts.append((t,min(all_throughputs)))
#     max_pts.append((t,max(all_throughputs)))
# var_pts = np.array(var_pts)
# min_pts = np.array(min_pts)
# max_pts = np.array(max_pts)


# p = Plot("","Number of threads", "Throughput")
# p.add("Var", *(var_pts.T))
# p.add("Min", *(min_pts.T))
# p.add("Max", *(max_pts.T))
# p.plot()
# exit()

