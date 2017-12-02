# Global variables
import sys, os, time, pickle, webbrowser
import numpy as np
# Get relative path to the data and this directory
home = os.path.expanduser("~")
print("Using home:", home)
sys.path += [home+"/Dropbox/Sync/bin"]
print("Using path:", sys.path)
from plotly_interface import *


sys.path.append(home + "/Git_Analytics/Analytics-Vis/Code")
curr_dir = os.path.abspath(os.path.dirname(__file__))
from settings import data_for_spec, spec_str_to_dict

NUM_PLOT_POINTS = 1000
NUM_SUBSAMPLES = 1000
NUM_TRIALS = 2
NUM_MODES = 3

sys.path += [home + "/Git_Analytics/Analytics-Vis/Data/VirginiaTech/Many_Trials/Data/"]
from plot_throughput import read_raw_file

# Code for generating the view file for a requested view
if __name__ == "__main__":
    # If there is an existing report, store it as 'old'
    os.system("mv report.html old-report.html")

    # Get the three primary command line arguments
    group_id = "VirginiaTech"
    project_names = ["Raw_CFQ_NOOP_64_Fread", "F_R_Size_Raw", "Many_Trials"]
    data_req_base = "[F_size={0}][R_Size={1}][Hyp_Sched=CFQ][VM_Sched=NOOP][Freq=2700000][Mode=Fread][Throughput=-inf_&le;_inf]"

    # Ordered in terms of largest to smallest discrepancy variance
    f_r_sizes = [(64,32),(1024,32),
                         (1024,512)]

    study = "Throughput"
    yaxis = "Percent Less"
    names = ["First",
             "Second",
             "True",
             # "All Together"
    ]

    data_file_names = [home+"/Git_Analytics/Analytics-Vis/Data/VirginiaTech/Raw_CFQ_NOOP_64_Fread/Data/Clean_Data.txt",
                       home+"/Git_Analytics/Analytics-Vis/Data/VirginiaTech/F_R_Size_Raw/Data/Clean_Data.txt",
                       home+"/Git_Analytics/Analytics-Vis/Data/VirginiaTech/Many_Trials/Data/Clean_Data.txt"
    ]

    # Cycle through the file and record sizes, generating a plot for each
    for f_size, r_size in f_r_sizes:
        # Get all three data sets from file
        data_req = data_req_base.format(f_size, r_size)
        old_out = sys.stdout; sys.stdout = open("/dev/null",'w')
        concat_data = []
        all_data = []
        for project in project_names:
            data, _,_,_ = data_for_spec(
                group_id, project, spec_str_to_dict(data_req))
            # We only actually want the study data from each of these
            all_data.append(data[study])
            concat_data += list(all_data[-1])
            print("="*70)
        sys.stdout = old_out
        concat_data = np.array(concat_data)

        # Generate important relavant meta-data about each set of points
        cdf_diffs = [cdf_second_diff(all_data[0]),
                     cdf_second_diff(all_data[1]),
                     cdf_second_diff(all_data[2])
        ]
        min_max = [(min(all_data[0]), max(all_data[0])),
                   (min(all_data[1]), max(all_data[1])),
                   (min(all_data[2]), max(all_data[2])),
        ]        
        true_min_max = (min(map(min,min_max)), max(map(max,min_max)))
        diff_min = min(min(diff[:,1]) for diff in cdf_diffs)
        diff_max = max(max(diff[:,1]) for diff in cdf_diffs)


        # Create all of the plots
        for f_name,n,mm,d,second_diff in zip(
                data_file_names,names,min_max,all_data,cdf_diffs):

            color_ind = names.index(n)
            n = "'%s'"%n
            # Create some identifying names based on current settings
            plot_name = "Fsize-%i_Rsize-%i"%(f_size, r_size)
            title = "Analysis of %s %s"%(n,plot_name)
            # Initialize the plot for displaying the CDF's
            pdf_plot = Plot(title, "Probabilty of Occurrence", study)
            # Produce plot of raw througput
            raw_plot = Plot(title.replace("PDF",study), "Trial Number", study)
            # Produce plot of the mode information
            mode_plot = Plot(title.replace("PDF","Mode"), "CDF''", study)

            color = pdf_plot.color( color_ind )

            #      Generate raw data plot     
            # ================================
            print("Generating raw data plot for %s..."%n)
            read_raw_file(raw_plot, f_name, f_size, r_size, n + " Raw Trials",
                          shade=False, color=color)

            #      Generate PDF Plot     
            # ===========================
            print("Generating 100 bin PDF for %s..."%n)
            bin_size = (mm[1]-mm[0])/100
            pdf_plot.add(n + " 100 Bin PDF", y_values=d, plot_type="histogram",
                         opacity=0.7, autobiny=False, fill_color=color,
                         histnorm='probability',
                         ybins=dict(start=mm[0],end=mm[1],size=bin_size))


            #      Generate Mode Plot     
            # ============================
            print("Generating mode bars for %s..."%n)
            modes = mode_points(second_diff)

            # Plot the top modes based on change in second derivative
            color = mode_plot.color( color_ind, brightness=1.0, alpha=0.8 )
            mode_x_vals = []
            mode_y_vals = []
            for i,pt in enumerate(modes[:NUM_MODES]):
                mode_x_vals += [0,     1,     None]
                mode_y_vals += [pt[0], pt[0], None]                
            mode_plot.add("%s Top %i Modes"%(n,NUM_MODES),
                          mode_x_vals, mode_y_vals, mode="lines",
                          color=color)

            # Plot the remaining modes
            color = mode_plot.color( color_ind, brightness=1.0, alpha=0.1 )
            mode_x_vals = []
            mode_y_vals = []
            for i,pt in enumerate(modes[NUM_MODES:]):
                mode_x_vals += [0,     1,     None]
                mode_y_vals += [pt[0], pt[0], None]                
            mode_plot.add("%s Remaining %i Modes"%
                          (n,len(modes)-NUM_MODES), mode_x_vals,
                          mode_y_vals, mode="lines",
                          color=color)

            # Plot the second difference plot
            m_color = mode_plot.color( color_ind, alpha=0.8 )
            modes = [list(pt) for pt in modes]
            modes.sort(key=lambda pt: pt[0])
            modes = np.array(modes)
            mode_plot.add("Normalized Absolute CDF''", modes[:,1],
                          modes[:,0], mode='markers',
                          #fill='tozerox', fill_color=color,
                          line_width=1, color=m_color, marker_size=3,
                           shade=False, show_in_legend=False)


            # Get all of the figures and multiplot them!
            raw_fig = raw_plot.plot(show=False, html=False)
            pdf_fig = pdf_plot.plot(show=False, html=False,
                                    layout=dict(barmode='overlay'))
            mode_fig = mode_plot.plot(show=False, html=False)
            multiplot([mode_fig, pdf_fig, raw_fig],
                      shared_y=True, show=False)
            os.system("cat temp-plot.html >> report.html")
            with open("report.html",'a') as f:
                f.write('\n\n<p style="page-break-before: always;"></p>\n\n')
            time.sleep(0.2)
        
    
    print("Opening %s in web browser."%(curr_dir + "/report.html"))
    webbrowser.open(curr_dir + "/report.html")
