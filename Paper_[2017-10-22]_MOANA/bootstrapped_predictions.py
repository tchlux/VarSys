import pickle, os
import numpy as np
from scipy.interpolate import splrep, splev
from util.algorithms import Delaunay
from util.data import Data
from util.plot import Plot, multiplot

# Path to raw text data
HOME_DIR = os.path.expanduser("~")
# DATA_PATH = os.path.join(HOME_DIR, *["Git","Git_Analytics", "Analytics-Vis",
#                                      "Data", "VirginiaTech",
#                                      "Many_Trials", "Data",
#                                      "Clean_Data.txt"])

DATA_PATH = os.path.join(HOME_DIR, *["Git", "Git_Analytics", "Analytics-Vis",
                                     "Data", "VirginiaTech",
                                     "Raw_Sampled", "Full_Data",
                                     "Clean_Data.txt"])

# Paths to local cache files (for speed of re-execution)
LOCAL_FILE = "IOzone_data_2016.pkl"
ORGANIZED_DATA_FILE = "data_dict.pkl"
# LOCAL_FILE = "Many_Trials_IOzone_data_2016.pkl"
# ORGANIZED_DATA_FILE = "many_trials_data_dict.pkl"
PREDICTED_RUNS_FILE = "predicted_runs.pkl"
NUM_RUNS_DICT = "%i_runs_%i_convergence_(%s,%s,%s).pkl"

# Known hard-coded data settings
OS_SCHEDULER = "Hyp_Sched"
VM_SCHEDULER = "VM_Sched"
FILE_SIZE = "F_size"
RECORD_SIZE = "R_Size"
THREAD_COUNT = "Threads"
TEST_TYPE = "Mode"
CPU_FREQUENCY = "Freq"
IO_THROUGHPUT = "Throughput"
# Data columns broken up by type of interest
CATEGORICAL_COLS = [OS_SCHEDULER, VM_SCHEDULER, TEST_TYPE]
CONTINUOUS_COLS = [FILE_SIZE, RECORD_SIZE, THREAD_COUNT, CPU_FREQUENCY]
TARGET_COLUMN = [IO_THROUGHPUT]
# Print settings (for showing more or less specific information)
SHOW_VALUES = False
PLOT_BOOTSTRAPPED = False
PLOT_CONVERGENCE = True
PLOT_CONVERGENCE_EXAMPLE = False
COLLECT_RUNS_DATA = False
PREDICT_RUNS = COLLECT_RUNS_DATA and True
RUN_PREDICTION_TRAIN_PERCENT = .9
RUN_PREDICTION_TEST_PERCENT = 1 - RUN_PREDICTION_TRAIN_PERCENT
RUN_PREDICTION_TRIALS = 1000

CONVERGENCE_PERCENT = .9

np.random.seed(2)

# Given a list of (i samples, lower bound, upper bound) that have
# different starting and ending widths with (near) monotone
# convergence from one to the other, return the number of "i samples"
# such that the (convergence / time) is maximized (most optimally
# spent time).
def calc_num_runs(i_low_upp):
    starting_width = i_low_upp[0][2] - i_low_upp[0][1]
    ending_width = i_low_upp[-1][2] - i_low_upp[-1][1]
    convergence = []
    time_convergence = []
    for i in range(len(i_low_upp)):
        width = i_low_upp[i][2] - i_low_upp[i][1]
        amount_converged = abs(width - starting_width) / abs(ending_width - starting_width)
        convergence.append( amount_converged )
        time_convergence.append( amount_converged / i_low_upp[i][0] )
    return convergence, time_convergence

# Function for calculating two bootstrapped standard deviation percentiles
def bootstrapped_std(values, sample_size=40, num_iters=2000,
                     lower_percentile=5, upper_percentile=95):
    stdevs = []
    for i in range(num_iters):
        stdevs.append( np.std(np.random.choice(values, sample_size)) )
    return ( np.percentile(stdevs,lower_percentile),
             np.percentile(stdevs, upper_percentile) )

# ==========================
#      Data preparation     
# ==========================

if not os.path.exists(LOCAL_FILE):
    print("Reading in the data...")
    data = Data.load(DATA_PATH, verbose=True).to_struct()
    with open(LOCAL_FILE,"wb") as f:
        print("Writing data to local pickle file..")
        pickle.dump(data,f)
else:
    print("Reading pickled data...")
    with open(LOCAL_FILE,"rb") as f:
        data = pickle.load(f)

print("Done loading data.")
if SHOW_VALUES:
    for h in list(data.dtype.names):
        if "Throughput" in h: continue
        print(h," "*(15 - len(h)),list(np.unique(data[h])))
    print()

if not os.path.exists(ORGANIZED_DATA_FILE):
    print("Getting unique continuous coordinates for modelling...")
    model_data = np.array([list(r) for r in
                           np.unique(data[CONTINUOUS_COLS].copy())])
    print("Moving data to dictionary...")
    lookup_data = {tuple(row):[] for row in np.unique(data[CATEGORICAL_COLS+CONTINUOUS_COLS])}
    print("Creating lists of throughputs associated with each config...")
    for row in data[CATEGORICAL_COLS+CONTINUOUS_COLS+TARGET_COLUMN]:
        row = list(row)
        config = tuple(row[:-1])
        target = row[-1]
        lookup_data[ config  ] += [target]
    print("Saving to file...")
    with open(ORGANIZED_DATA_FILE, "wb") as f:
        pickle.dump((model_data, lookup_data), f)
else:
    print("Reading model and lookup data from file...")
    with open(ORGANIZED_DATA_FILE, "rb") as f:
        (model_data, lookup_data) = pickle.load(f)

print("Done.")

# ===========================
#      Data manipulation     
# ===========================

# Using the default files, collect the data on the number of runs
# needed at each system configuration
if COLLECT_RUNS_DATA:
    os_sched = "CFQ"
    vm_sched = "NOOP"
    test_type = "Fread"
    print(model_data, len(model_data))
    if not os.path.exists(PREDICTED_RUNS_FILE):
        # Holders for the values
        i_low_upp = {}
        num_runs_needed = {}
        for row in model_data:
            config = (os_sched, vm_sched, test_type) + tuple(row)
            target = lookup_data[config]
            i_low_upp[config] = []
            for i in range(10, len(target)+1, 1):
                print("\rBootstrapping %i"%i, end="")
                low,upp = bootstrapped_std(target, sample_size=i)
                i_low_upp[config].append( (i,low,upp) )
            convergence, convergence_by_time = calc_num_runs(i_low_upp[config])
            num_runs_needed[config] = (convergence, convergence_by_time)
            print()

        # Save data to file
        with open(PREDICTED_RUNS_FILE, "wb") as f:
            pickle.dump((i_low_upp, num_runs_needed), f)
    else:
        with open(PREDICTED_RUNS_FILE, "rb") as f:
            i_low_upp, num_runs_needed = pickle.load(f)
    
    num_runs_dict = NUM_RUNS_DICT%(list(i_low_upp.values())[0][-1][0],
                                   int(CONVERGENCE_PERCENT*100),
                                   os_sched, vm_sched, test_type)

    if not os.path.exists(num_runs_dict):
        num_runs = {}
        for key in num_runs_needed:
            convergence, convergence_by_time = np.asarray( num_runs_needed[key] )
            iters = np.array([i for (i,low,upp) in i_low_upp[key]])
            conv_index = [i for i,v in enumerate(convergence) if v >= CONVERGENCE_PERCENT][0]
            if PLOT_CONVERGENCE_EXAMPLE:
                print(key)
                p = Plot("Convergence of confidence interval by number of runs",
                         "Number of runs",
                         "Ratio of estimated interval width to final interval width")

                p.add("90% Convergence", [iters[0], iters[-1]], [.9,.9], mode='lines')
                # func = lambda x: (np.sin(x/10)/50) + .9
                # p.add_func("",func, [iters[0], iters[-1]])
                p.add("Convergence", iters, convergence)
                annotation_x = iters[conv_index]
                annotation_y = convergence[conv_index]
                p.add_annotation("90% Convergence", annotation_x, annotation_y,
                                 ax=30,ay=20,y_anchor='top')
                p.plot(show_legend=False, file_name="NEW Convergence_Stoping_Point_%i.html"%(iters[-1]))
                exit()
            else:
                num_runs[key] = conv_index

        print("Saving '%s'..."%(num_runs_dict))
        with open(num_runs_dict, "wb") as f:
            pickle.dump(num_runs, f)
    else:
        print("Loading '%s'..."%(num_runs_dict))
        with open(num_runs_dict, "rb") as f:
            num_runs = pickle.load(f)


# Build a model that predicts the number of runs needed, compare the
# model's performance with leave-some-out.
if (COLLECT_RUNS_DATA and PREDICT_RUNS):
    print("Predicting runs")
    runs_guess_file = num_runs_dict.replace("convergence_", "convergence_guesses_")
    if not os.path.exists(runs_guess_file):
        run_prediction_guesses = {}
        # Generate the target values for the predictive model
        target = np.asarray([
            num_runs[(os_sched, vm_sched, test_type) + tuple(row)]
                for row in model_data], dtype=float)
        # Run a series of prediction trials that have randomly selected 
        for i in range(RUN_PREDICTION_TRIALS):
            print("\rPredicting trial %i of %i"%(i, RUN_PREDICTION_TRIALS), end="")
            train_size = int(round(len(model_data) * RUN_PREDICTION_TRAIN_PERCENT))
            test_size = len(model_data) - train_size
            indices = list(range(len(model_data)))
            np.random.shuffle(indices)
            # Separate model data into training and testing sets randomly
            train_data = np.asarray(model_data[indices[:train_size]], dtype=float)
            test_data = np.asarray(model_data[indices[train_size:]], dtype=float)
            # Normalize the training data in order to make dimensiosn
            # equally weighted when predicting
            shift = np.min(train_data, axis=0)
            train_data -= shift
            test_data -= shift
            scale = np.max(train_data, axis=0)
            train_data /= scale
            test_data /= scale
            # Now make the predictions
            train_target = target[indices[:train_size]]
            test_target = target[indices[train_size:]]
            model = Delaunay()
            model.fit(train_data, train_target)
            guesses = model(test_data)
            # Now store the estimates for each of the testing points
            for row,estimate in zip(test_data, guesses):
                row = (row * scale) + shift
                config = (os_sched, vm_sched, test_type) + tuple(int(round(v)) for v in row)
                if config not in run_prediction_guesses:
                    run_prediction_guesses[config] = []
                run_prediction_guesses[config].append( estimate )


        # Save the results
        print("Saving '%s'..."%(runs_guess_file))
        with open(runs_guess_file, "wb") as f:
            pickle.dump(run_prediction_guesses, f)
    else:
        print("Loading '%s'..."%(runs_guess_file))
        with open(runs_guess_file, "rb") as f:
            run_prediction_guesses = pickle.load(f)

    predicted = []
    actual = []
    for key in run_prediction_guesses:
        guesses = [round(v) for v in run_prediction_guesses[key]]
        predicted += guesses
        actual += [num_runs[key]] * len(guesses)
    predicted = np.asarray(predicted)
    actual = np.asarray(actual)
    p = Plot("Predicted Runs Needed for Convergence versus Truth",
             "Predicted Number of Runs Needed",
             "True Number of Runs Needed")
    p.add("",predicted, actual, color=p.color(1))
    m,b = np.polyfit(predicted, actual, 1)
    func = lambda x: m*x + b
    p.add_func("", func, [min(predicted), max(predicted)], color=p.color(0))
    # p.plot(show_legend=False)
    errors = predicted - actual
    print("Most under:  ", int(round(min(errors))))
    print("Most over:   ", int(round(max(errors))))
    print("Median error:", int(round(np.median(errors))))
    print("Mean error:  ", int(round(np.mean(errors))))
    print("Error stdev: ", int(round(np.std(errors))))


if PLOT_CONVERGENCE:
    os_sched = "CFQ"
    vm_sched = "NOOP"
    test_type = "Fread"
    file_size = [256] # [64, 256, 1024]
    record_size = [128] # [32, 128, 512]
    threads = [16] # [1, 2, 4, 8, 16, 32, 64, 128, 256]
    frequencies = [2700000] # [1200000, 1400000, 1500000, 1600000, 1800000, 1900000, 2000000, 2100000, 2300000, 2400000, 2500000, 2700000, 2800000, 2900000, 3000000, 3001000]
    # Holders for the plotted values
    x_values = []
    lower = []
    upper = []
    guessed_lower = []
    guessed_upper = []
    # Tracking for automatically creating the right plot
    ordered_names = ["File Size", "Record Size", "Thread Count", "Frequency"]
    ordered_lists = [file_size, record_size, threads, frequencies]
    interest_index = np.argmax(list(map(len,ordered_lists)))
    print()
    print("Model Data")
    print(model_data)
    print()
    for row in model_data:
        if not ((row[0] in file_size) and
                (row[1] in record_size) and
                (row[2] in threads) and
                (row[3] in frequencies)):
            continue
        print("Found the right row")
        print()
        # Produce new data set without this point, generate an
        # estimated target distribution in order to bootstrap a guess
        new_data = np.array([tuple(r) for r in model_data if tuple(r) != tuple(row)])
        # Normalize the data
        shift = np.min(new_data, axis=0)
        new_data = new_data - shift
        scale = np.max(new_data, axis=0)
        new_data = new_data / scale
        # Build a model to identify the simplex points
        model = Delaunay()
        model.fit(new_data)
        points, weights = model._predict(((row-shift)/scale)[None,:])
        # Transform the points and weights into a matrix and single vector
        points = np.array([new_data[i] for i in points[0]])
        weights = weights[0]
        # De-normalize the points and extract only if non-zero weight
        points = np.array([[int(round(v)) for v in p*scale+shift]
                           for (p,w) in zip(points,weights) if w > 0])
        weights = np.array([w for w in weights if w > 0])
        # Generate estimated target values for this configuration
        guessed_target = None
        for p,w in zip(points, weights):
            config = (os_sched, vm_sched, test_type) + tuple(p)
            if type(guessed_target) == type(None):
                guessed_target = np.array(lookup_data[config]) * w
            else:
                guessed_target = guessed_target + np.array(lookup_data[config]) * w
        # Generate the lower / upper bounds of the true and guessed intervals
        config = (os_sched, vm_sched, test_type) + tuple(row)
        target = lookup_data[config]
        for i in range(10, len(target)+1, 1):
            print("\rBootstrapping %i"%i, end="")
            # Get the true values
            low,upp = bootstrapped_std(target, sample_size=i)
            x_values.append(i)
            lower.append(low)
            upper.append(upp)
            # Get the guessed values
            low,upp = bootstrapped_std(guessed_target, sample_size=i)
            guessed_lower.append(low)
            guessed_upper.append(upp)
        print()

    # =======================================
    #      Generate the convergence plot     
    # =======================================

    # Create plot of true convergence
    p = Plot("", "Sample Size", "I/O Stdev", font_family="times", font_size=20)
    guessed_p = Plot("90% Confidence Interval with Increasing Sample Size", 
                     "", "Predicted I/O Stdev", font_family="times", font_size=20)

    # Calculate some useful values
    min_max = [min(x_values), max(x_values)]
    color = p.color(1)
    annotation_color = p.color(0)
    annotation_font_size = 17

    # Add the convergence series
    p.add("", x_values, lower, color=color, mode="markers+lines",
          fill="tonexty", fill_opacity=0.2)
    p.add("", x_values, upper, color=color, mode="markers+lines")

    # Add the estimated convergence series
    guessed_p.add("", x_values, guessed_lower, color=color, mode="lines+markers", 
                  fill="tonexty", fill_opacity=0.2)
    guessed_p.add("", x_values, guessed_upper, color=color, mode="lines+markers")

    # =====================================
    #      Add convergence annotations     
    # =====================================

    for conv_perc in [0.5, 0.75, 0.9]:
        # Extract the convergence information
        convergence, convergence_by_time = calc_num_runs(list(zip(x_values,lower,upper)))
        conv_index = [i for i,v in enumerate(convergence) if v >= conv_perc][0]
        annotation_x = x_values[conv_index]
        annotation_y = lower[conv_index]

        ax = annotation_x + 1
        ay = annotation_y - 200000

        # Add vertical line for convergence
        p.add("",[annotation_x]*2, [lower[conv_index], upper[conv_index]],
              color=annotation_color, mode="lines")
        p.add_annotation("%.0f%% Convergence<br>to observed truth"%(100*(conv_perc)),
                         annotation_x, annotation_y, 
                         ax=ax, ay=ay, y_anchor="top", align="right",
                         font_size=annotation_font_size,
                         font_family="times", x_anchor="left")


        # Extract the guessed convergence information
        guessed_convergence, guessed_convergence_by_time = calc_num_runs(
            list(zip(x_values,guessed_lower,guessed_upper)))
        guessed_conv_index = [i for i,v in enumerate(guessed_convergence) 
                              if v >= conv_perc][0]
        guessed_annotation_x = x_values[guessed_conv_index]
        guessed_annotation_y = guessed_lower[guessed_conv_index]

        ax = guessed_annotation_x + .5
        ay = 14

        # Add vertical line for convergence
        guessed_p.add("",[guessed_annotation_x]*2, 
                      [guessed_lower[guessed_conv_index], 
                       guessed_upper[guessed_conv_index]],
                      color=annotation_color, mode="lines")
        guessed_p.add_annotation("%.0f%% Convergence Guess<br>Error of %.1f%% (%i)"%(
            100*(conv_perc),
            100*(guessed_annotation_x - annotation_x)/x_values[-1], 
            guessed_annotation_x - annotation_x),
                                 guessed_annotation_x, guessed_annotation_y, 
                                 ax=ax, ay=ay, y_anchor="top", align="right",
                                 font_size=annotation_font_size,
                                 font_family="times", x_anchor="left")

    # Create plot
    fig1 = p.plot(show_legend=False,
                  axis_settings=dict(showgrid=False), 
                  html=False)
    # Create plot
    fig2 = guessed_p.plot(show_legend=False,
                          axis_settings=dict(showgrid=False), 
                          html=False)


    multiplot([[fig1],[fig2]], shared_x=True, file_name="NEW double_convergence.html")



if PLOT_BOOTSTRAPPED:
    os_sched = "CFQ"
    vm_sched = "NOOP"
    test_type = "Fread"
    file_size = [64, 256, 1024]
    record_size = [32] # [32, 128, 512]
    threads = [16] # [1, 2, 4, 8, 16, 32, 64, 128, 256]
    frequencies = [2400000] # [1200000, 1400000, 1500000, 1600000, 1800000, 1900000, 2000000, 2100000, 2300000, 2400000, 2500000, 2700000, 2800000, 2900000, 3000000, 3001000]
    # Holders for the plotted values
    x_values = []
    lower = []
    upper = []
    # Tracking for automatically creating the right plot
    ordered_names = ["File Size", "Record Size", "Thread Count", "Frequency"]
    ordered_lists = [file_size, record_size, threads, frequencies]
    interest_index = np.argmax(list(map(len,ordered_lists)))
    for row in model_data:
        if not ((row[0] in file_size) and
                (row[1] in record_size) and
                (row[2] in threads) and
                (row[3] in frequencies)):
            continue
        config = (os_sched, vm_sched, test_type) + tuple(row)
        target = lookup_data[config]
        low,upp = bootstrapped_std(target)
        x_values.append(row[interest_index])
        lower.append(low)
        upper.append(upp)
    # Create nice plot
    p = Plot("NEW 90%% Confidence Interval of Standard Deviation versus %s"
             %(ordered_names[interest_index]),
             ordered_names[interest_index], 
             "Standard Deviation of I/O Throughput")
    min_max = [min(x_values), max(x_values)]
    print("Using splines that are order %i."%(min(3,len(x_values)-1)))
    fit = splrep(x_values, lower, k=min(3,len(x_values)-1))
    fit = lambda x_val, fit=fit: splev(x_val, fit)
    p.add("", x_values, lower, color=color)
    p.add_func("5th Percentile", fit, min_max, mode="lines",
               fill="tonexty", fill_opacity=0.2, color=color)
    fit = splrep(x_values, upper, k=min(3,len(x_values)-1))
    fit = lambda x_val, fit=fit: splev(x_val, fit)
    p.add_func("95th Percentile", fit, min_max, mode="lines",
               color=color)
    p.add("", x_values, upper, color=color)
    # Generate a description automatically
    desc = file_size+record_size+threads+frequencies
    for i in range(len(ordered_lists[interest_index])):
        desc.pop(interest_index)
    desc = ",".join(list(map(str,desc)))
    # Create plot
    p.plot(file_name="NEW_SD_90_vs_%s_(%s).html"%(
        ordered_names[interest_index].replace(" ","_"), desc),
           show_legend=False, axis_settings=dict(showgrid=False),)
           # layout=dict(plot_bgcolor="#fcfcfc",
           #             margin=dict(l=50,r=50,b=100,t=100)))



# print(len(model_data))
# print(len(lookup_data))
# print(max(len(v) for v in lookup_data.values()))
# print(min(len(v) for v in lookup_data.values()))



# print("Creating model...")
# model = Delaunay()
# model.fit(model_data)
# print(len(model_data))
# print(model_data[:10])
