import os
import numpy as np
from util.data import Struct
from util.plotly import Plot, multiplot
from util.stats import cdf_fit_func

test = "readers"

# Shade and color series
# 
# 0 -- "Frequency"      "initial_writers"
# 1 -- "File Size"      "rewriters"      
# 2 -- "Record Size"    "readers"        
# 3 -- "Num Threads"    "re-readers"     
#                       "random_readers" 
#                       "random_writers" 

# color_by = 1
# shade_by = 2
# reduction = 2800000
# reduction_by = "Frequency"

color_by = 0
shade_by = 3
reduction = 1048576 # 4
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


# Reduce the set of distributions based on file size
# dists = dists[ dists["File Size"] == file_size ]
dists = dists[ dists[reduction_by] == reduction ]

# # Reduce the number of threads visible
# dists = dists[ dists["Num Threads"] < 48 ]

# Summarize the data that was loaded...
dists[dists.names[:-1]].summarize(13)

unique_vals = dists.unique()
# Generate the lookup for defining colors
color_options = sorted(unique_vals[dists.names[color_by]])
colors = dict(zip(color_options, list(range(len(color_options)))))

# Generate the lookup for defining shades
shade_options = sorted(unique_vals[dists.names[shade_by]])
shades = dict(zip(shade_options, np.linspace(.5,1,len(shade_options))))

print("Generating plot...")
p = Plot(f"CDFs for {reduction_by} {reduction} (Colored by {dists.names[color_by]}, Shaded by {dists.names[shade_by]})", 
         "Log Throughput", "CDF Value", font_family="times")
# Add legend entries for each series
for c in color_options:
    name = f"{dists.names[color_by]} = {c}"
    c = colors[c]
    p.add(name, [None],[None], color=p.color(c), group=c, mode="lines")

# Add all of the CDFs to the plot
for row in dists:
    # Get the CDF for this row
    f = row[-1]
    # Get the color and shade
    c = colors[row[color_by]]
    b = shades[row[shade_by]]
    # Add plot series for this CDF
    p.add_func("", f, f(), color=p.color(c, brightness=b), group=c,
               show_in_legend=False, hoverinfo="text",
               text=str(row[:-1]))


# Put the legend on the bottom.
legend_settings = dict(
    xanchor = "center",
    yanchor = "top",
    x = .5,
    y = -.15,
    orientation = "h",
)
# Generate plot file
p.show(file_name=f"[2018-04-14]_({reduction_by.replace(' ','_')}_{reduction})_CDFs_({test}).html",
       fixed=False, x_axis_settings=dict(type="log"), legend=legend_settings)
