# Global variables
import sys, os, time, pickle, webbrowser
import numpy as np
# Get relative path to the data and this directory
home = os.path.expanduser("~")
print("Using home:", home)
sys.path += [home+"/Dropbox/Sync/bin"]
print("Using path:", sys.path)
from plotly_interface import *


sys.path = [home + "/Git_Analytics/Analytics-Vis/Code"] + sys.path
curr_dir = os.path.abspath(os.path.dirname(__file__))
from settings import data_for_spec, spec_str_to_dict

NUM_PLOT_POINTS = 1000
NUM_SUBSAMPLES = 1000
NUM_TRIALS = 2

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

        # Create some identifying names based on current settings
        plot_name = "Fsize-%i_Rsize-%i"%(f_size, r_size)
        title = "CDF Comparison of %s"%(plot_name)
        # Initialize the plot for displaying the CDF's
        cdf_plot = Plot(title, study, yaxis)
        pdf_plot = Plot(title.replace("CDF","PDF"),
                        "Probabilty of Occurrence", study)

        # Produce plot of raw througput
        raw_plot = Plot(title.replace("CDF",study), "Trial Number", study)
        for f_name,n in zip(data_file_names,names):
            read_raw_file(raw_plot, f_name, f_size, r_size, n,
                          shade=False, group=n)
            print(n)
        raw_fig = raw_plot.plot(show=False)
        # We will combine this with a vertical histogram soon

        # Generate the initial functions for the two old samples
        funcs = [cdf_fit_func(all_data[0]),
                 cdf_fit_func(all_data[1]),
                 cdf_fit_func(all_data[2]),
                 # cdf_fit_func(concat_data)
        ]
        min_max = [(min(all_data[0]), max(all_data[0])),
                   (min(all_data[1]), max(all_data[1])),
                   (min(all_data[2]), max(all_data[2])),
                   # (min(concat_data), max(concat_data))
        ]
        
        # Plot the plain CDF's
        true_min_max = (min(map(min,min_max)), max(map(max,min_max)))
        for f,mm,n,d in zip(funcs, min_max, names, all_data):
            line_color = None if names != names[2] else 'rgb(10,10,10)'
            line_width = None if n != names[2] else 1.0
            cdf_plot.add_func(n, f, true_min_max,
                              line_color=line_color,
                              line_width=line_width)

            # Simultaneously build the PDF plot for the raw data
            bin_size = (mm[1]-mm[0])/100
            pdf_plot.add(n, y_values=d, plot_type="histogram",
                         opacity=0.7, autobiny=False,
                         histnorm='probability',group=n,show_in_legend=False,
                         ybins=dict(start=mm[0],end=mm[1],size=bin_size))


        pdf_fig = pdf_plot.plot(show=False,html=False,
                                layout=dict(barmode='overlay'))

        multiplot([pdf_fig, raw_fig], [[0.0,0.15],[0.2,1.0]],
                  shared_y=True, show=False)

        os.system("cat temp-plot.html >> report.html")
        time.sleep(0.2)
        
        # Generate the data for the percentiles (by subsampling a lot)
        subsamples = NUM_SUBSAMPLES
        subsample_size = 40
        percentiles = [0,10,20,30,40,60,70,80,90,100]
        # Get the min and max of the master set
        x_points = np.linspace(min_max[2][0], min_max[2][1], NUM_PLOT_POINTS)
        cdf_sampled_funcs = [
            cdf_fit_func(np.random.choice(all_data[2], size=subsample_size))
            for i in range(subsamples)]
        y_points = np.array([fit(x_points) for fit in cdf_sampled_funcs])
        # Generate extra functions for the percentiles of the large sample
        perc_names = ["%ith Percentile"%p for p in percentiles]
        # Generate the functions for each of the percentiles of the large sample
        perc_funcs = percentile_funcs(x_points, y_points, percentiles, min_max[2])

        plot_percentiles(cdf_plot, names[2], perc_funcs, percentiles,
                         true_min_max)

        cdf_plot.plot(show=False)
        os.system("cat temp-plot.html >> report.html")
        time.sleep(0.2)


        #      Compare the KS-Differences for the CDF's     
        # ==================================================

        # all_true_ks_diff = ks_diff(cdf_fit_func(concat_data), funcs[2])
        one_two_ks_diff = ks_diff(funcs[0], funcs[1])
        one_true_ks_diff = ks_diff(funcs[0], funcs[2])
        two_true_ks_diff = ks_diff(funcs[1], funcs[2])

        confidence_alpha = lambda alpha: np.sqrt(-np.log(alpha/2)/2)
        ks_pass = lambda ks_diff, n1, n2, alpha=0.01: ks_diff > (
            confidence_alpha(alpha) * np.sqrt((n1 + n2) / (n1*n2)))
        smallest_alpha = lambda diff, n1, n2: 2 * np.exp(
            (diff / np.sqrt((n1 + n2) / (n1*n2)))**2 * -2)

        print("Confidence(alpha) with alpha = 0.01:", confidence_alpha(0.01))
        print()

        print("Sample sizes:")
        # print("","All    ",len(concat_data))
        print("","First: ",len(all_data[0]))
        print("","Second:",len(all_data[1]))
        print("","True:  ",len(all_data[2]))
        print()

        print("KS-Differences:")
        # print("","All 500 -- True 420:  ", all_true_ks_diff)
        print("","First 40 -- Second 40:", one_two_ks_diff)
        print("","First  40 -- True 420:", one_true_ks_diff)
        print("","Second 40 -- True 420:", two_true_ks_diff)
        print()

        print("KS-Differences Pass at alpha 0.01 (significantly different):")
        # print("","All 500 -- True 420:  ",ks_pass(
        #     all_true_ks_diff,len(concat_data), len(all_data[2])))
        print("","First 40 -- Second 40:",ks_pass(
            one_two_ks_diff,len(all_data[0]), len(all_data[1])))
        print("","First  40 -- True 420:",ks_pass(
            one_true_ks_diff,len(all_data[0]), len(all_data[2])))
        print("","Second 40 -- True 420:",ks_pass(
            two_true_ks_diff,len(all_data[1]), len(all_data[2])))
        print()

        print("Largest possible confidence in KS-Difference:")
        # print("","All 500 -- True 420:  ", 100.0*(1.0 - smallest_alpha(            
        #     all_true_ks_diff,len(concat_data), len(all_data[2]))))
        print("","First 40 -- Second 40:", 100.0*(1.0 - smallest_alpha(
            one_two_ks_diff,len(all_data[0]), len(all_data[1]))))
        print("","First  40 -- True 420:", 100.0*(1.0 - smallest_alpha(
            one_true_ks_diff,len(all_data[0]), len(all_data[2]))))
        print("","Second 40 -- True 420:", 100.0*(1.0 - smallest_alpha(
            two_true_ks_diff,len(all_data[1]), len(all_data[2]))))
        print()

        print("Smallest possible alpha for passing KS-Difference:")
        # print("","All 500 -- True 420:  ",smallest_alpha(
        #     all_true_ks_diff,len(concat_data), len(all_data[2])))
        print("","First 40 -- Second 40:", smallest_alpha(
            one_two_ks_diff,len(all_data[0]), len(all_data[1])))
        print("","First  40 -- True 420:", smallest_alpha(
            one_true_ks_diff,len(all_data[0]), len(all_data[2])))
        print("","Second 40 -- True 420:", smallest_alpha(
            two_true_ks_diff,len(all_data[1]), len(all_data[2])))
        print()


        #      Test the convergence of the CDF's     
        # ===========================================

        print("Testing Convergence..")
        trials = []
        subsample_sizes = np.array(list(range(2, len(all_data[2])+1 )))
        # for test in range(NUM_TRIALS):
        #     print(" testing",test)
        #     one_trial = []
        #     for size in subsample_sizes:
        #         fit = cdf_fit_func(np.random.choice(all_data[2], size=size))
        #         one_trial.append(ks_diff(fit, funcs[2]))
        #     trials.append(one_trial)

        # percentiles = [0,10,20,30,40,50,60,70,80,90,100]
        # trials = np.array(trials)

        # # Save the trials to file (for later analysis
        # with open("Trials_"+title+".pkl", "wb") as f:
        #     pickle.dump(trials, f)

        with open("First_100_Trial_Report/Trials_"+title+".pkl", "rb") as f:
            trials = pickle.load(f)

        min_max = (min(subsample_sizes), max(subsample_sizes))
        convergence_funcs = percentile_funcs(subsample_sizes, trials,
                                             percentiles, min_max)
        perc_names = ["%ith Percentile"%p for p in percentiles]
        center_color = np.array([50,100,255,0.7])
        outer_color = np.array([150,200,255,0.4])
        
        conv_plot = Plot("Convergence of " + title, "Num Samples",
                          "Max CDF Diff (KS Score)")

        # Plot the percentiles of the convergence
        plot_percentiles(conv_plot, "", convergence_funcs, percentiles,
                         min_max, center_color, outer_color)
        
        conv_plot.plot(show=False)

        # Produce a temporary plot file (do not open automatically)
        os.system("cat temp-plot.html >> report.html")
        time.sleep(0.2)        

        alpha_plot = Plot("Alpha with increasing samples for " + title,
                          "Num Samples", "Alpha (likelihood of coming from same distribution)")

        trials = [[smallest_alpha(row[i], size, len(all_data[2]))
                   for i,size in enumerate(subsample_sizes)]
                  for row in trials]
        trials = np.array(trials)
        convergence_funcs = percentile_funcs(subsample_sizes, trials,
                                             percentiles, min_max)

        plot_percentiles(alpha_plot, "", convergence_funcs,
                         percentiles, min_max, center_color,
                         outer_color)
        alpha_plot.plot(show_legend=False, show=False)        
        # Produce a temporary plot file (do not open automatically)
        os.system("cat temp-plot.html >> report.html")
        time.sleep(0.2)

        # title = "<strong>File Size %i — Record Size %i</strong>"%(f_size, r_size)
        # print(REPORT_ENTRY%(title, one_ks_diff, two_ks_diff, file_name))
            
    
    webbrowser.open("report.html")

# Run command:
# python3 code.py > report.html && google-chrome report.html