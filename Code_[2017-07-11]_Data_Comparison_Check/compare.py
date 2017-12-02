import pickle, sys, os, time, webbrowser
import numpy as np
from Resources.plotly_interface import *

CURR_DIR = os.path.abspath(os.path.curdir)
FILE_DELIM = ","
RESOURCE_DIR = os.path.join(CURR_DIR, "Resources")
OLD_DATA = os.path.join(RESOURCE_DIR,"VarSys-2016.txt")
OLD_DATA_NP = os.path.join(RESOURCE_DIR, "VarSys-2016.pkl")
LOAD_RAW_OLD_DATA = False # False -> use saved pkl file if it exists
DATA_TYPE_KEYS = ["Machine", "Store", "Journal", "Hyp", "Hyp_Sched",
                  "VM_Sched", "RCU", "F_Size", "R_Size", "Threads",
                  "Mode", "Freq", "Throughput"]
MEAN_OPACITY = 0.3
DATA_TYPES = {
    "Machine":'S3',
    "Store":'S3',
    "Journal":'S3',
    "Hyp":'S3',
    "Hyp_Sched":'S4',
    "VM_Sched":'S4',
    "RCU":float,
    "F_Size":float,
    "R_Size":float,
    "Threads":float,
    "Mode":'S13',
    "Freq":float,
    "Throughput":float
}

# Read the help string for the program
with open("README.txt") as f: HELP_STRING = f.read()

HTML_INFO_ALERT = '''
<div class="alert alert-info alert-dismissable">
  <a href="#" class="close" data-dismiss="alert" aria-label="close">&times;</a>
  <strong>%s</strong>%s
</div>
'''

# Takes: .format(title, mean_diff, stdev_diff, var_diff, fig_str)
HTML_BEGIN = '''
<!DOCTYPE html>
<html lang="en">
  <head>
    <title>Report on F Size and R Size Variance</title>
    <!-- Default setup for a webpage with basic bootstrap -->
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <link rel="stylesheet" href="{0}">
    <link rel="stylesheet" href="{1}">

    <script src="{2}"></script>
    <script src="{3}"></script>

    <style type="text/css">
        #Brand {{ cursor: default; }}
        #Project {{ float: right; }}
        #Project-Name:hover {{ color: #555; }}
        #Project-Name {{
            float: right;
            font-size: 16px;
            padding: 14px 5px;
            color: #777;
            text-decoration: none;
            cursor: default;
        }}
        div.entry {{      
            border-top: 1px solid #eee;
            border-bottom: 1px solid #eee;
            background-color: #fdfdfd;
            width: 90vw;
            margin: auto;
            margin-top: 10px;
            margin-bottom: 10px;
        }}
        div.entry-title {{
            padding: 5px 0px;
            font-size: 14px;
        }}
        div.holder {{
            padding: 10px;
            height: 45vw;
        }}
        div.entry-text {{
            font-size: 14px;
            float: left;
            display: inline-block;
            padding-top: 15vw;
            width: 25vw;
            height: 45vw;
        }}
        figure {{
            height: 100%;
            width: 60vw;
            float: right;
            display: inline-block;
        }}
        figure.plain {{
            height: (45vw - 20px);
            width: 100%
        }}
    </style>
  </head>

  <body id="Body">
    <!-- Nav bar -->
    <div class="nav">
	<a id="Brand">Virginia Tech</a>
	<div id="Project">
	  <a id="Project-Name">VarSys</a>
	</div>
    </div>

'''.format(os.path.join(RESOURCE_DIR,"bootstrap.min.css"),
           os.path.join(RESOURCE_DIR,"analysis.css"),
           os.path.join(RESOURCE_DIR,"jquery.min.js"),
           os.path.join(RESOURCE_DIR,"bootstrap.min.js"))

REPORT_ENTRY = '''
    <div class="entry">
      <div class="entry-title" onclick="$(this.parentNode.children[1]).collapse('toggle')">
	%s
      </div>
      <div class="holder collapse in" aria-expanded="true">
	<div class="entry-text">
	  <p style="text-decoration: underline; font-size: 16px; color: #333;">Statistical differences</p>
	  <strong>Mean:</strong> %s
	  <br>
	  <strong>Stdev:</strong> %s
	  <br>
	  <strong>Variance:</strong> %s
	</div>
	<figure>
	  %s
	</figure>
        <br style="clear: both">
      </div>      
    </div>
'''

PLAIN_ENTRY = '''
    <div class="entry">
      <div class="entry-title" onclick="$(this.parentNode.children[1]).collapse('toggle')">
	%s
      </div>
      <div class="holder collapse in" aria-expanded="true">
        %s
        <br style="clear: both">
      </div>      
    </div>
'''

# ;setTimeout(function(){window.dispatchEvent(new Event('resize'))},500);
# onload="$(this.parentNode.children[1]).collapse('toggle')"

PLAIN_ENTRY_COLLAPSED = '''
    <div class="entry">
      <div class="entry-title" onclick="{$(this.parentNode.children[1]).collapse('toggle');setTimeout(function(){window.dispatchEvent(new Event('resize'))},500);}">
	%s
      </div>
      <div class="holder collapse" aria-expanded="false">
        %s
        <br style="clear: both">
      </div>      
    </div>
'''

HTML_END = '''
  </body>
</html>
'''

# Return the data types of the columns based on the headers of a data file
def get_dtypes(f_name):
    dtypes = []
    with open(f_name) as f:
        header = [h.strip() for h in f.readline().split(FILE_DELIM)]
        for h in header:
            if h not in DATA_TYPES:
                print()
                print("ERROR: Header '%s' does not exist in original data."%h)
                exit()
            else:
                dtypes.append( (h,DATA_TYPES[h]) )
    return dtypes

# =============================
#      Main Execution Code     
# =============================

if __name__ == "__main__":
    # Check for the correct number of command line arguments
    if len(sys.argv) != 2:
        print(HELP_STRING)
        exit()

    existing_file = os.path.exists(OLD_DATA_NP)
    if LOAD_RAW_OLD_DATA or (not existing_file):
        old_dtypes = get_dtypes(OLD_DATA)
        print("NOTICE: Loading old data set from file for comparison...")
        if not LOAD_RAW_OLD_DATA:
            print("        This action takes some time (~60 seconds), and should only")
            print("        occur once. It generates the file 'Resources/VarSys-2016.pkl'")
        old_data = np.loadtxt(OLD_DATA, dtype=old_dtypes,
                              delimiter=FILE_DELIM, skiprows=1)
        print("        Done loading old data set.")

        print("NOTICE: Dumping old data into pickle file for faster loading..")
        with open(OLD_DATA_NP, "wb") as f:
            pickle.dump(old_data, f)
        print("        Done dumping old data.")

    print("NOTICE: Loading old data from file...")
    with open(OLD_DATA_NP, "rb") as f:
        old_data =pickle.load(f)
    print("        Done loading old data.")

    new_data_filename = sys.argv[1]
    print("NOTICE: Loading data from file '%s'..."%new_data_filename)
    new_dtypes = get_dtypes(new_data_filename)
    new_data = np.loadtxt(new_data_filename, dtype=new_dtypes,
                          delimiter=FILE_DELIM, skiprows=1)
    # If the new data only had 1 row, reshape it to be a matrix
    if len(new_data.shape) == 0:
        new_data = new_data.reshape((1,))
    print("        Done loading '%s'."%new_data_filename)

    # Make the pretty version of each file name
    new_data_name = os.path.basename(new_data_filename)
    new_data_name = new_data_name[:-1-new_data_name[::-1].index(".")]
    old_data_name = os.path.basename(OLD_DATA)
    old_data_name = old_data_name[:-1-old_data_name[::-1].index(".")]

    # Check to make sure the new data has throughput
    if "Throughput" not in new_data.dtype.names:
        print()
        print("ERROR:  No column labeled 'Throughput' in new data.")
        print()
        exit()

    print("NOTICE: Identifying comparable rows from old data...")
    # Get the controllable settings in the old and new data
    all_settings = list(old_data.dtype.names)
    all_settings.remove("Throughput")
    settings = list(new_data.dtype.names)
    settings.remove("Throughput")

    # Next, identify the remaining rows that have the same settings
    to_keep = np.array([False]*old_data.shape[0])
    for i,row in enumerate(new_data):
        matches_row = np.array([True]*old_data.shape[0])
        for s in settings:
            matches_row = np.logical_and(matches_row, old_data[s] == row[s])
            # Error checking for values that don't exist in the old data
            if not np.any(matches_row):
                print(" "*30)
                print("ERROR:  No matching historical data for settings")
                max_len_name = max(map(len,settings))
                for h in settings:
                    val = row[h]
                    if (type(val) == np.bytes_): val = str(val)[1:]
                    spaces = max_len_name - len(h) + 3
                    print("         %s%s%s"%(h," "*spaces,val))
                print()
                exit()
        to_keep = np.logical_or(to_keep, matches_row)
        print("        %0.2f%%"%(100.0*(i+1)/new_data.shape[0]),end="\r")
    print("        Reducing old data...")
    old_data = old_data[to_keep]
    print("        Done reducing old data, found %s points."%len(old_data))

    print("NOTICE: Creating plots..")
    # - old data as series (shaded regions for mean, 3 stdevs)
    #   new data spread across same width (only add the mean line)
    # - 100 bin PDF old-new comparison plot, 2 series each partially transparent
    # - For each header, do colored series for the options (same
    #   colors different symbols for two data sets)
    
    report_name = "%s_report[%s].html"%(new_data_name,time.strftime("%Y-%m-%d"))
    with open(report_name,"w") as f:
        print(HTML_BEGIN, file=f)

        # Warnings related to data reduction
        for s in all_settings:
            if s not in settings:
                unique_vals = np.unique(old_data[s])
                if len(unique_vals) > 1:
                    print(HTML_INFO_ALERT%("WARNING!",
                        (" No settings for '%s' specified in new data. "+
                         "Displaying old data for all %i known options.")%(s,len(unique_vals)))
                          ,file=f)


        #      First report on the raw data     
        # ======================================
        title = "<strong>%s</strong>"%("Raw Data")
        # Get the statistical differences between the old and new data
        old, new = old_data["Throughput"].var(), new_data["Throughput"].var()
        var_diff = 100.0 * abs(old - new) / old
        old, new = old_data["Throughput"].std(), new_data["Throughput"].std()
        std_diff = 100.0 * abs(old - new) / old
        old, new = old_data["Throughput"].mean(), new_data["Throughput"].mean()
        mean_diff = 100.0 * abs(old - new) / old
        # Convert the statistical differences to strings
        mean_diff = "%.2f%%"%mean_diff
        std_diff = "%.2f%%"%std_diff
        var_diff = "%.2f%%"%var_diff
        if len(new_data) == 1: std_diff = var_diff = 'N/A'
        # Plot the differences
        p = Plot("","Collection Order","Throughput")
        old_x = np.array(list(range(len(old_data)))) + 1.0
        new_x = np.array(list(range(len(new_data)))) + 1.0
        # Add the new data and its mean line
        p.add(new_data_name, new_x, new_data["Throughput"], group=new_data_name, shade=False)
        color = p.color(p.color_num, alpha=MEAN_OPACITY)
        p.add(new_data_name+" Mean", [1,max(old_x[-1], new_x[-1])],[new,new],
              group=new_data_name, color=color, mode="lines", show_in_legend=False)
        # Add the old data and its mean line
        p.add(old_data_name, old_x, old_data["Throughput"], group=old_data_name, shade=False)
        color = p.color(p.color_num, alpha=MEAN_OPACITY)
        p.add(old_data_name+" Mean", [1,max(old_x[-1], new_x[-1])],[old,old],
              group=old_data_name, color=color, mode="lines", show_in_legend=False)
        # Add the upper and lower bounds of the old data as dashed lines
        old_upper = max(old_data["Throughput"])
        p.add(old_data_name+" std-up", [1,max(old_x[-1], new_x[-1])],[old_upper, old_upper],
              group=old_data_name, color=color, dash=10, mode="lines", show_in_legend=False)
        old_lower = min(old_data["Throughput"])
        p.add(old_data_name+" std-dn", [1,max(old_x[-1], new_x[-1])],[old_lower, old_lower],
              group=old_data_name, color=color, dash=10, mode="lines", show_in_legend=False)
        # Generate the plot and add the output to the report file
        p.plot(show=False)
        with open("temp-plot.html") as i: file_str = i.read()
        print(REPORT_ENTRY%(title, mean_diff, std_diff, var_diff, file_str),file=f)


        #      Second report on the 100-bin PDFs     
        # ===========================================
        title = "<strong>%s</strong>"%("Probability Distributions")
        p = Plot("","Throughput Value","Probability of Occurrence")
        min_val, max_val = min(old_data["Throughput"]), max(old_data["Throughput"])
        bin_size = (max_val - min_val) / 100
        p.add(new_data_name, x_values=new_data["Throughput"], 
              plot_type="histogram", opacity=0.7, autobinx=False,
              histnorm='probability',
              xbins=dict(start=min_val,end=max_val,size=bin_size))
        p.add(old_data_name, x_values=old_data["Throughput"], 
              plot_type="histogram", opacity=0.7, autobinx=False,
              histnorm='probability',
              xbins=dict(start=min_val,end=max_val,size=bin_size))
        # Generate the plot and add the output to the report file
        p.plot(show=False)
        with open("temp-plot.html") as i: file_str = i.read()
        print(PLAIN_ENTRY%(title, file_str),file=f)


        #      Third (et al) report(s) on the different column values     
        # ================================================================
        for s in all_settings:
            title = "<strong>%s Analysis</strong>"%(s.replace("_"," "))
            # Make a plot comparing the new and old versions of this setting
            p = Plot("","Collection Order","Throughput")
            p.color_num += 1
            # Add the new data series only if there are multiple options
            options = np.unique(old_data[s])
            if len(options) == 1: continue
            for o in options:
                name = new_data_name + " " + str(o)
                # Get only the relevant data
                if s in list(new_data.dtype.names):
                    y_vals = new_data[np.where(new_data[s]==o)]["Throughput"]
                else:
                    y_vals = []
                x_vals = list(range(len(y_vals)))
                p.add(name, x_vals, y_vals, group=o, shade=False,
                      symbol="square", show_in_legend=len(y_vals) > 0)
                # Add the old data (different symbol)
                color = p.color(p.color_num)
                name = old_data_name + " " + str(o)
                y_vals = old_data[np.where(old_data[s]==o)]["Throughput"]
                x_vals = list(range(len(y_vals)))
                p.add(name, x_vals, y_vals, group=o, color=color, shade=False, symbol="circle")
            # Generate the plot and add the output to the report file
            p.plot(show=False)#,layout={"width":"90vw", "height":"45vw"})
            with open("temp-plot.html") as i: file_str = i.read()
            print(PLAIN_ENTRY_COLLAPSED%(title, file_str),file=f)

        print(HTML_END,file=f)

    print("        Done, opening report saved as")
    print("         '%s'."%report_name)
    print()

    # Remove the temporary plot in this folder
    os.remove(os.path.join(CURR_DIR, "temp-plot.html"))
    path_to_file = os.path.join(os.path.abspath(os.path.curdir),report_name)

    print("IGNORE: Any remaining notifications are produced by python webbrowser.")
    print()
    webbrowser.open("file://"+path_to_file)



