from util.data import Data
from util.plot import Plot
import pygal, os

data_path = os.path.join(os.path.expanduser("~"),"Downloads",
                         "Stats_Data_October_2018.csv")
data = Data.load(data_path)

# ['ID', 'Media', 'Frequency', 'File.Size', 'Record.Size',
#  'Num.Threads', 'Test', 'runs', 'distcomp', 'disname']
all_names = ['Frequency', 'File.Size', 'Record.Size', 
             'Num.Threads', 'Test', 'disname']
x_name = "Num.Threads"
normalized = False
plot_type = "normalized" if normalized else "count"
output_file = f"[2018-10-09]_Distnum_Visual_({plot_type}).html"

# Initialize the output file to be empty.
with open(output_file, "w") as f:  f.write("")


# Make all of the plots.
for x_name in all_names:
    y_name = "distcomp"
    # Get the unique x-values
    x = sorted(set(data[x_name]))
    y = sorted(set(data[y_name]))
    print(x_name,"--",x)
    # Get the total counts of data instances with each x value.
    total_counts = {v:sum(1 for i in data[x_name]==v) for v in x}
    # Initialize the plot
    plot = pygal.StackedBar()
    plot.title = f"Observed Number of Distributions by {x_name}"
    plot.x_title = x_name
    plot.y_title = "Percentage" if normalized else "Count"
    plot.x_labels = list(map(str, x))

    # Get the unique y-values associated with each x and plot them
    for dist_count in y:
        sub_data = data[data[y_name] == dist_count]
        vals = []
        for v in x:
            sub_sub_data = sub_data[sub_data[x_name] == v]
            # Compute the stacked value
            if not normalized: val = len(sub_sub_data)
            else:              val = 100. * len(sub_sub_data) / total_counts[v]
            # Add the appropriate value to the stack.
            if len(sub_sub_data) == 0: vals.append( None )
            else:                      vals.append( val )
        # plot.add(f"{y_name} -- {dist_count}", vals)
        plot.add(f"{dist_count} Mode{'s' if dist_count > 1 else ''}", vals)

    string = plot.render()
    with open(output_file, "a") as f:
        f.write(str(string, "utf-8"))

print("Opening plot in web browser..")
import webbrowser
path_name = "file://" + os.path.join(os.path.abspath("."), output_file)
webbrowser.open(path_name)

