import numpy as np
import pickle, sys, pathlib, os
from util.multi_dim_analysis import make_test_data, MDA_Iterator
from notifications import send_email
from util.decorators import timeout_after
import util.algorithms
from py_settings import *

# 1) using all 4 continuous dimensions (file size, record size,
#    threads, frequency), normalized to be inside the unit
#    (hyper)cube, predict variance and then mean for a single
#    selection of categorical settings (CFQ hypervisor scheduler, NOOP
#    virtual machine Scheduler, File-read test). Potentially modify
#    these categorical settings later.

# 1a) perform leave-1-out testing for each datum *inside the convex
#     hull* of the 4D data 

# 1b) perform a sliding range of (95-5, 90-10, 85-15, ...) 
#     training-testing splits of the data ensuring *all* testing data
#     is inside the convex hull of the training data. Using N trials
#     for each split, where N is determined by the smallest number of
#     trials that gives us 95% confidence in the predicted results for
#     *all* splits (assuming normally distributed results) 

# 1c) train and test each algorithm with each of the splits from above
#     (each algorithm is given the same input and asked to predict the
#     same values)

# 2) repeat step 1 using all 4 pairs of 3 continuous dimensions,
#    averaging the variance / mean for overlapping configurations 

# 3) repeat step 1 using all 6 pairs of 2 continuous dimensions,
#    averaging the variance / mean for overlapping configurations 

# 4) repeat step 1 using all 4 continuous dimensions, averaging the
#    variance / mean for overlapping configurations 

SOURCE_DIR = "Data_Source"

DATA = os.path.join(SOURCE_DIR, "CFQ_NOOP_Fread.pkl")
TARGET = "Throughput"
# Two data sets, one for mean throughput data, the other for variance
# of throughput data.
MEAN_MDA_DIR = "VarSys-2016-Mean_MDA_Data_0"
VAR_MDA_DIR = "VarSys-2016-Variance_MDA_Data_0"

# DATA = os.path.join(SOURCE_DIR, "Stage_1_Fread.pkl")
# TARGET = "Throughput"
# # Two data sets, one for mean throughput data, the other for variance
# # of throughput data.
# MEAN_MDA_DIR = "VarSys-2017-Mean_MDA_Data_0"
# VAR_MDA_DIR = "VarSys-2017-Variance_MDA_Data_0"

# Load the original data file (and get the column names)
with open(DATA,"rb") as f:
    data = pickle.load(f)
    settings = [v for v in data.dtype.names if len(np.unique(data[v])) > 1]
    settings.remove(TARGET)

# Generate the MDA data set if necessary
if ((not os.path.exists(MEAN_MDA_DIR)) or
    (not os.path.exists(VAR_MDA_DIR))):
    # Filter down the data into the unique x-coordinates
    unique_pts = {}
    for row in data:
        # Get the unique system settings and the target value
        x_point = tuple(row[s] for s in settings)
        val = row[TARGET]
        # If the unique point already exists, simply append to the
        # list of recorded values, if it does not, initialize a list
        if x_point not in unique_pts:
            unique_pts[x_point] = [val]
        else:
            unique_pts[x_point].append(val)
    # Convert the unique points into a list and sort them
    data_keys = list(unique_pts.keys())
    data_keys.sort()

    # Generate the mean and variance data by concatenating the statistic
    # value of all data points associated with each unique key
    stat_data = []
    for key in data_keys:
        vals = np.array(unique_pts[key])
        stat_data.append( key + (vals.mean(), vals.var()) )
    # List of the column names and their data types, simply adding two
    # columns for mean and variance for conversion back into numpy
    # structured arrays.
    dtypes = [(s,float) for s in settings] + [("Mean",float), ("Variance", float)]
    stat_data = np.array(stat_data, dtype=dtypes)

    print()
    print("Controllable settings:\n  %s"%(settings))
    print("  with %s raw data points"%len(data))
    print("  with %s unique data points"%len(stat_data))
    print()

    # Generate the testing data files if necessary
    inputs = stat_data[settings]
    if not(os.path.exists(MEAN_MDA_DIR)):
        output = stat_data["Mean"]
        make_test_data(inputs, output, "VarSys-2016"+"-Mean", settings)
    if not(os.path.exists(VAR_MDA_DIR)):
        output = stat_data["Variance"]
        make_test_data(inputs, output, "VarSys-2016"+"-Variance", settings)

    print("Ready to test and evaluate")


# Generate the list of algorithms (from all possible algorithms in util)
algorithms = []
for name in dir(util.algorithms):
    try:
        if issubclass(getattr(util.algorithms,name), util.algorithms.Approximator):
            algorithms.append(getattr(util.algorithms,name)())
    except:
        pass

# All algorithms have:
#  <alg>.fit(<inputs>, <responses>)
#  <alg>(<inputs>) -> <responses>

# The names of the outputs from testing
data_header = settings.copy() + ["Algorithm", "Dimension",
                                 "Number_of_Trials",
                                 "Train_Percentage",
                                 "Number_of_Points",
                                 "Number_of_Guesses", 
                                 "Mean_Fit_Time",
                                 "Mean_Evaluation_Time", 
                                 "Truth",
                                 "Mean_Guess", 
                                 "Guess_Min",
                                 "Guess_25th",
                                 "Guess_50th",
                                 "Guess_75th", 
                                 "Guess_Max", 
                                 "Mean_Error", 
                                 "Absolute_Mean_Error",
                                 "Relative_Mean_Error"] 

# String to recommend certain usage
HELP_STRING = """
USAGE:
  python3 %s <test_type> <alg1> [<alg2>, ..., <algn>]

Where <test_type> is one of "mean" or "variance", and each <alg> is
one of the algorithms provided in the algorithms module of the util
package. A full set of tests will be run on the appropriate MDA data
set.
"""%(__file__)

if len(sys.argv) < 3:
    print(HELP_STRING)
    exit()

to_predict = (sys.argv[1]).title()
to_test = sys.argv[2:]
for alg_name in to_test:
    print(alg_name)
    # Retreive the requested algorithm
    for alg in algorithms:
        name = GET_CLASS_NAME(alg)
        if (name == alg_name): break
    else:
        print("ERROR:   Couldn't find '%s' in known algorithms."%alg_name)
        continue

    # Notify user, create output file name
    print()
    print("Testing: '%s', "%(name), end="")
    output_data_folder = "Output-"+name+"_"+to_predict+"-"+OUTPUT_DATA_FILE[:-4]
    # Create an output data folder
    os.makedirs(output_data_folder)
    output_data_file = "Output-"+name+"_"+to_predict+"-"+OUTPUT_DATA_FILE
    checkpoint_file = "Checkpoint-%s.txt"%(name)
    print("results being stored in '%s'"%output_data_file)
    print()

    # Initialize the results file if it does not exist
    if (not os.path.exists(output_data_file)):
        with AtomicOpen(output_data_file, "w") as f:
            print(",".join(data_header), file=f)
    # Initialize a checkpoint file if it does not exist
    if (not os.path.exists(checkpoint_file)):
        pathlib.Path(checkpoint_file).touch()
        checkpoints = ""
    else:
        with AtomicOpen(checkpoint_file) as f:
            checkpoints = f.read()
    
    # Get a file iterator that will cycle through all of the
    # train/test files for evaluation. Initialize variables.
    if "variance" in sys.argv:
        file_iterator = MDA_Iterator(VAR_MDA_DIR)
    else:
        file_iterator = MDA_Iterator(MEAN_MDA_DIR)
    data = {}
    start = time.time()
    # "estimated time remaining" using linear extrapolation
    etr = lambda p: 100.0 * ((time.time() - start) / p)
    # Initialize holders for writing data to file
    alg_performance_data = {}
    number_of_trials = 0
    # Cycle through the train/test files in the Multi Dimensional Analysis
    for (train, test, dim, cols, num) in file_iterator:
        # Skip files that have already been checkpointed
        if (os.path.dirname(train)) in checkpoints: continue

        # Writing results to file is done after all of the train-test
        # files have been iterated from one ratio of training and
        # testing data volume. A new set is encountered when num == 0.
        if (num == 0) and (number_of_trials > 0): 
            lines_of_info = []
            # Collect and average all of the testing data over the
            # unique test points (from multiple trials)
            unique_vals = sorted(alg_performance_data.keys())
            # Collect some default info that is the same for each unique point
            common_info = ([name, str(input_dimension)] +
                           list(map( CLEAN_NUMBER_STRING,
                [number_of_trials, train_perc, train_size+test_size])))

            # Save the raw data into a directory heirarchy
            split_folder = os.path.basename(os.path.dirname(train))
            col_folder = os.path.basename(os.path.dirname(os.path.dirname(train)))
            raw_folder = os.path.join(output_data_folder, col_folder, split_folder)
            os.makedirs(raw_folder)
            for x_point in unique_vals:
                raw_file_name = "-".join(x_point)
                guesses = np.array(alg_performance_data[x_point]["Guess"])
                with open(os.path.join(raw_folder, raw_file_name), "w") as f:
                    for g in guesses: print(g, file=f)

            # Cycle through each of the unique 'true' data points,
            # collect statistics on that data point, record to file
            for x_point in unique_vals:
                # Initialize information on this unique data point
                info = list(x_point) + common_info
                # Get the lists of information of relevant algorithm performance
                guesses = np.array(alg_performance_data[x_point]["Guess"])
                num_guesses = len(guesses)
                fit_times = np.array(alg_performance_data[x_point]["Fit_Time"])
                eval_times = np.array(alg_performance_data[x_point]["Evaluation_Time"])
                truth = alg_performance_data[x_point]["Truth"]
                guess_min = np.min(guesses)
                guess_max = np.max(guesses)
                mean_error = (guesses - truth).mean()
                mean_abs_error = (abs(guesses - truth)).mean()
                rel_error = (((guesses - truth) / truth) if truth != 0 
                             else ((guesses - truth) / SMALL_NUMBER)).mean()
                # Condense all the numeric statistics into the
                # shortest length strings possible
                info += list(map( CLEAN_NUMBER_STRING,
                    [num_guesses, fit_times.mean(), eval_times.mean(),
                     truth, guesses.mean(), guess_min,
                     np.percentile(guesses,25),
                     np.percentile(guesses,50), 
                     np.percentile(guesses,75), guess_max,
                     mean_error, mean_abs_error, rel_error] ))
                # add this line of performance data as a string
                lines_of_info.append( ",".join(info) )

            # Write the information to file
            with AtomicOpen(output_data_file, "a") as f:
                for line in lines_of_info:
                    print(line, file=f)
            # Checkpoint the completion of this training file
            with AtomicOpen(checkpoint_file, "a") as f:
                print(working_dir, file=f)
            # Reset the stored performance data for next set of tests
            alg_performance_data = {}

        # Testing an algorithm on a given train / test pair
        # - load the training and testing data files
        # - train the algorithm on the training data
        # - test the algorithm on the testing data
        perc_complete = 100.0*file_iterator.current/len(file_iterator)

        print("%.2f%% of %i tests"%(perc_complete, len(file_iterator)),
              "(%.0fs -> %.0fs)"%(time.time()-start,etr(perc_complete)),
              end="\r", flush=True)

        train_data = np.loadtxt(train, dtype=float, delimiter=",", ndmin=2)
        test_data = np.loadtxt(test, dtype=float, delimiter=",", ndmin=2)
        # Extract the points and values from the data file
        train_points = train_data[:,:-1]
        train_values = train_data[:,-1]
        test_points = test_data[:,:-1]
        test_values = test_data[:,-1]
        # Remove these variables from memory (they are not needed anymore)
        del(train_data, test_data)
        # Normalize the training and testing data
        data_shift = np.min(train_points, axis=0)
        train_points -= data_shift
        test_points -= data_shift
        data_scale = np.max(train_points, axis=0)
        train_points /= data_scale
        test_points /= data_scale
        # try:
        # Fit the training points
        # Get any extra arguments to pass to individual alrogithm fits
        extra_args = EXTRA_ARGS.get(name,[])
        # Wrap the fit function to monitor timeouts
        fit_func = timeout_after(FIT_TIMEOUT_SECONDS,FAIL_VALUE)(alg.fit)
        # Time the fit operation
        fit_time = time.time()
        fit_output = fit_func(train_points, train_values, *extra_args)
        fit_time = time.time() - fit_time
        if ((type(fit_output) == type(FAIL_VALUE)) and
            (fit_output == FAIL_VALUE)):
            print()
            print('ERROR: %s failed to fit within %i seconds.'%(
                name,FIT_TIMEOUT_SECONDS))
            print("Training file: '%s'"%(train))
            print("Testing file:  '%s'"%(test))
            print()
            break # Cancel the rest of data collection for this algorithm
        # Wrap the approximation function to monitor timeouts
        approx_func = timeout_after(APPROX_TIMEOUT_SECONDS,FAIL_VALUE)(alg)
        # Evaluate at the testing points
        eval_time = time.time()
        test_approx = approx_func(test_points)
        eval_time = time.time() - eval_time
        if ((type(test_approx) == type(FAIL_VALUE)) and
            (test_approx == FAIL_VALUE)):
            print()
            print('ERROR: %s failed to create approximation within %i seconds.'%(
                name, APPROX_TIMEOUT_SECONDS))
            print("Training file: '%s'"%(train))
            print("Testing file:  '%s'"%(test))
            print()
            break # Cancel the rest of data collection for this algorithm
        # except:
        #     print()
        #     print('ERROR: %s failed to fit and approximate data.'%(name))
        #     print("Training file: '%s'"%(train))
        #     print("Testing file:  '%s'"%(test))
        #     print()
        #     break # Cancel the rest of data collection for this algorithm            

        # De-normalize the data so that correct values are writtne to file
        train_points *= data_scale
        test_points *= data_scale
        train_points += data_shift
        test_points += data_shift

        # If this is the first test in a new round, initialize
        # statistics on the size of training and testing sets
        if (len(alg_performance_data) == 0):
            total_pts = len(train_points) + len(test_points)
            train_size = len(train_points)
            test_size = len(test_points)
            train_perc = "%.2f"%(100.0*train_size/total_pts)
            test_perc = "%.2f"%(100.0*test_size/total_pts)

        # Collect the lines of information from this test to record
        for i in range(len(test_points)):
            # Retreive the 'i'th test point as a data point and get
            # the approximation produced by the algorithm
            x = test_points[i]
            if len(test_points) > 1:
                approx = test_approx[i]
            else:
                approx = test_approx
            # Collect all relevant performance information on this
            # test point in this test train/test pair
            row = []
            x_index = 0
            # Cycle the controllable settings, (creates unique points)
            for col_name in settings:
                # See if that setting is currently being tested
                # 'cols' is the list of settings being tested
                if (col_name in cols):
                    # Collect all settings from the data into 'row'
                    row.append( CLEAN_NUMBER_STRING(x[x_index]) )
                    x_index += 1
                else:
                    # Since this setting is not being testing, append '-1'
                    row.append( "-1" )
            # Convert row into an immutable form of data
            row = tuple(row)
            # Initialize a holder for this test point if necessary
            if (row not in alg_performance_data):
                truth = test_values[i]
                alg_performance_data[row] = {
                    "Fit_Time":[], "Evaluation_Time":[],
                    "Guess":[], "Truth":truth }
            # Append test-run information for this test point
            alg_performance_data[row]["Guess"].append(approx)
            alg_performance_data[row]["Fit_Time"].append(fit_time)
            alg_performance_data[row]["Evaluation_Time"].append(eval_time)

        # Continue updating batch size, the last update will be the
        # true size of the batch of train/test files
        number_of_trials = num + 1
        # Store the working directory of this train / test pair
        working_dir = os.path.dirname(train)
        # Store the dimension (in case of save)
        input_dimension = dim
    print()


if len(to_test) > 0:
    # Notify me that the process has been completed
    send_email(subject="Algorithm evaluation complete")

# find <dir> -name "<name>" -exec rm -rf "{}" \;
