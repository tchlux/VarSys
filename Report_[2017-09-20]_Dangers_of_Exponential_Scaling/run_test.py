import pickle, os
from util.data import read_struct
import numpy as np
# Add a field to the unique columns that will hold the throughput std
from numpy.lib.recfunctions import append_fields

data_file = "Li_Results.txt"
train_pkl = "Train.pkl"
test_pkl = "Test.pkl"
train = "Train.txt"
test = "Test.txt"
target = "Throughput"
numerical_cols = ['F_size', 'R_Size']
categorical_cols = ['Hyp_Sched', 'VM_Sched', 'Mode']
testing_f_r_sizes = [(512,32),(512,128),(512,256),(768,32),(768,128)]
output_file = "New_Predictions.csv"



# SUMMARY:
#    Return a structured array with only the unique combinations of
#    desired columns and a summary statistic describing the set of
#    values that all share similar settings.
# 
# INPUTS:
#    data: 
#       A sorted numpy structured array (lexigraphic by column)
#    non_target_cols:
#       A list of the column headers that should be used to identify 
#       unique rows
#    target:
#       A target column header for a numerical column in "data".
#    statistic: (Default = np.std)
#       A function that converts a numpy array of numbers into a
#       single number. 
#    verbose: (Default = True)
#       A boolean indicating whether or not print statements about
#       progress should be produced during the conversion.
# 
# RETURNS:
#    A sorted numpy structured array (same order as "data")  with the
#    "target" column representing a summary statistic of the repeated
#    settings and an additional column, 'Count' that gives the number
#    of points encountered at each individual unique setting. 
def reduce_to_unique_vals(data, non_target_cols, target,
                          statistic=np.std, verbose=True):
    unique_data = np.unique(data[non_target_cols])
    unique_data = append_fields(unique_data, target,
                                np.zeros(unique_data.shape[0]))
    unique_data = append_fields(unique_data, "Count",
                                np.zeros(unique_data.shape[0], dtype=int))
    # Initialize indices for collecting sets of identical rows
    group_start_idx = 0
    group_end_idx = 0
    unique_current_idx = 0
    unique_setting = unique_data[non_target_cols][unique_current_idx]
    # Coalesce all data into the statistics
    while (group_end_idx < len(data)):
        if verbose: print("\r%i / %i"%(group_end_idx, len(data)), end="")
        if not (data[non_target_cols][group_end_idx] ==
                unique_data[non_target_cols][unique_current_idx]):
            # Find the standard deviation of this group of numbers
            unique_data[target][unique_current_idx] = (
                data[target][group_start_idx:group_end_idx].std())
            # Record the number of runs that were found for this setting
            unique_data["Count"][unique_current_idx] = (
                group_end_idx - group_start_idx)
            # Increment the current unique search pattern
            unique_current_idx += 1
            # Update the starting index of the current group, and the setting
            group_start_idx = group_end_idx
            unique_setting = unique_data[non_target_cols][unique_current_idx]
        # Increment the ending index of the group
        group_end_idx += 1
    if verbose: print()
    return unique_data


if __name__ == "__main__":
    print("Loading data...")

    # Load in the testing data
    if not (os.path.exists(test_pkl) and os.path.exists(train_pkl)):
        test_data = read_struct(test, verbose=True)
        # Reduce to only the columns with multiple unique values
        to_keep = []
        reduced = {}
        for header in test_data.dtype.names:
            unique_values = np.unique(test_data[header])
            if len(unique_values) > 1:
                to_keep.append( header )
            else:
                reduced[header] = unique_values[0]
        test_data = test_data[to_keep]

        # Identify the non-target columns (for unique identification)
        non_target_cols = [h for h in test_data.dtype.names if h != target]

        #      Convert the data into distribution information     
        # ========================================================
        # unique_test_data = reduce_to_unique_vals(
        #     test_data, non_target_cols, target, verbose=True)
        unique_test_data = np.unique(test_data[non_target_cols])
        unique_test_data = append_fields(unique_test_data, target,
                                    np.zeros(unique_test_data.shape[0]))
        for i in range(len(unique_test_data)):
            unique_test_data[target][i] = np.std(test_data[target][i*40:(i+1)*40])
        test_data = unique_test_data

        print()
        print("Reduced:",reduced)
        print("To Keep:",to_keep)
        print()

        train_data = read_struct(train, verbose=True)    
        print("Reducing training data..")
        # Make sure to grab same single settings from testing
        for key in reduced:
            train_data = train_data[train_data[key] == reduced[key]]
        # Reduce the training data to those columns relevant to testing
        train_data = train_data[to_keep]
        #      Convert the data into distribution information     
        # ========================================================
        # unique_train_data = reduce_to_unique_vals(
        #     train_data, non_target_cols, target, verbose=True)
        unique_train_data = np.unique(train_data[non_target_cols])
        unique_train_data = append_fields(unique_train_data, target,
                                    np.zeros(unique_train_data.shape[0]))
        for i in range(len(unique_train_data)):
            unique_train_data[target][i] = np.std(train_data[target][i*40:(i+1)*40])
        train_data = unique_train_data

        print("Saving pickle files..")
        # Save to pickle files
        with open(test_pkl, "wb") as f:
            pickle.dump(test_data, f)
        with open(train_pkl, "wb") as f:
            pickle.dump(train_data, f)
    else:
        with open(test_pkl, "rb") as f:
            test_data = pickle.load(f)
        with open(train_pkl, "rb") as f:
            train_data = pickle.load(f)
    print("","done")
    print()

    # # Reduce testing data to those specific rows
    # testing_indices = [tuple(row) in testing_f_r_sizes
    #                    for row in test_data[numerical_cols]]
    # # Extract the settings for known values from the test data
    # known_values = [tuple(r) for r in np.unique(train_data[numerical_cols])]
    # known_test_indices = [tuple(row) in known_values
    #                       for row in test_data[numerical_cols]]
    # known_test_data = test_data[known_test_indices]
    # # Reduce the test data to those interesting points
    # test_data = test_data[testing_indices]

    print()
    print("Headers:           ", list(test_data.dtype.names))
    print("Training size:     ",train_data.shape[0])
    print("Testing size:      ",test_data.shape[0])
    print()

    from util.algorithms import Delaunay
    from util.algorithms import LSHEP
    from util.algorithms import MARS
    from util.algorithms import MaxBoxMesh
    from util.algorithms import BayesTree
    from util.algorithms import CLASS_NAME

    algorithms = [Delaunay, LSHEP, MARS, MaxBoxMesh, BayesTree]
    error_data = [["Mode", "IO_Scheduler", "VM_Scheduler", "File_Size",
                    "Record_Size", "Original_Std", "True_Std"]]
    error_data[0] += list(map(lambda c: CLASS_NAME(c()) + "", algorithms))
    print(len(error_data), error_data)

    # Cycle unique settings and get the trianing / testing data
    i = 0
    unique_categoricals = np.unique(train_data[categorical_cols])
    for triple in unique_categoricals:
        i += 1
        print("\r%i / %i"%(i,len(unique_categoricals)),end="")
        # Extract the three categories from the triple
        io, vmio, mode = triple
        # Get the training data
        current_training_indices = train_data[categorical_cols] == triple
        current_training = train_data[current_training_indices]
        # Get the testing data
        current_testing_indices = test_data[categorical_cols] == triple
        current_testing = test_data[current_testing_indices]
        # Convert training and testing into standard numpy arrays
        train_x = np.array(list(map(
            list,current_training[numerical_cols])), dtype=float)
        train_y = np.array(current_training[target], dtype=float)
        test_x = np.array(list(map(
            list,current_testing[numerical_cols])), dtype=float)
        test_y = np.array(current_testing[target],  dtype=float)
        # Fit the predictive models
        predictors = [alg() for alg in algorithms]
        for p in predictors: p.fit(train_x, train_y)
        # Cycle through the testing points and add data to file
        for ((f_size,r_size), actual) in zip(test_x,test_y):
            # Generate predictions
            x_pt = np.array([f_size, r_size])
            pred_y = [alg(x_pt) for alg in predictors]
            # Identify the old value (if it exists)
            repeated_train = np.all([f_size, r_size] == train_x, axis=1)
            old_val = train_y[repeated_train][0] if any(repeated_train) else 0.0
            # Generate data row
            new_data_row = [mode, io, vmio, int(f_size), int(r_size),
                            old_val, actual] + pred_y
            error_data.append(new_data_row)
    print()

    print("Saving prediction results to file...")
    with open(output_file, "w") as f:
        for row in error_data:
            print(",".join(list(map(str,row))), file=f)

    error_data = read_struct(output_file)
    print(error_data)

