from util.data import Struct

data = Struct().load("prediction_results.csv")

data.summarize(15)
strange_data = data[ data[data.names[-1]] == 0 ]
print()
strange_data.summarize()
print(strange_data)

dists = Struct().load("data_readers_[cdfs].dill.gz")
dists = dists[ dists["File Size"] == 1048576 ]

print(dists)
dists[dists.names[:-1]].summarize()

print("Generating plot...")
from util.plotly import Plot
brights = {1200000:0.5,
           1600000:0.625,
           2000000:0.75,
           2300000:0.875,
           2800000:1.0
}
brights = {
    4    :.5,
    8    :.54,
    16   :.58,
    32   :.62,
    64   :.66,
    128  :.70,
    256  :.74,
    512  :.78,
    1024 :.82,
    2048 :.86,
    4096 :.90,
    8192 :.95,
    16384:1.0,
}

shown = set()
colors = []
p = Plot("CDFs for Large File Sizes (Colored by Thread, Brightness by Record Size)", "Throughput", "CDF Value")
for row in dists:
    f = row[-1]
    # Get threads (and color)
    threads = row[-2]
    name = f"{threads} Thread{'s' if (threads > 1) else ''}"
    if (threads not in colors): 
        colors.append(threads)

    c = colors.index(threads)
    # Get frequency (and brightness)
    bright = brights[row[2]]
    show_in_legend = (bright == 1.0) and (name not in shown)
    p.add_func(name, f, f(), color=p.color(c, brightness=bright),
               group=c, show_in_legend=show_in_legend,
               hoverinfo="text", text=name+" "+str(row[:-1]))
    if show_in_legend: shown.add(name)
p.show(file_name="[2018-04-14]_Large_File_Size_CDFs.html")
