from util.algorithms import Delaunay
from util.data import read_struct
import numpy as np

# SOURCE_COLUMN_NAMES = ["Machine","Media","Frequency","File Size",
#                        "Record Size","Num Threads","Iteration","Test",
#                        "Throughput","Runtime"]
# SOURCE_COLUMN_DTYPES = [str, str, int, int, int, int, int, str, float, int]

# CONFIG_COLUMN_NAMES = ["Frequency", "File Size", "Record Size", 
#                        "Num Threads", "Iterations"]
# CONFIG_COLUMN_DTYPES = [int, int, int, int, int]

MODELLED_COLUMNS = ["Frequency", "File Size", "Record Size", "Num Threads"]
TARGET_COLUMN = "Runtime"
ITERATION_COLUMN = "Iterations"

HELP_MESSAGE = """
Program for estimating the runtime of a given set of system
configurations. Expected usage:
  $ python3 estimate_runtime.py <config file name>.csv
  > X days, Y hours, and Z minutes.

The provided configuration file should be a comma separated values
file with 5 columns of data. The five columns in the file should
correspond with "Frequency" (HZ), "File Size" (KB), "Record Size" (KB), 
"Num Threads", and "Iterations" respectively. The file should look as follows:

<BEGIN FILE>
Frequency,File Size,Record Size,Num Threads,Iterations
3500000,12288,4,4,100
3500000,12288,4,16,100
2600000,12288,1024,64,200
1400000,24576,16384,64,50
... <more rows continued> ...
<END OF FILE>

The output of the program will be the estimated compute time of the
selected configurations based on performance data collected from
IoZone in Octoboer of 2017.
"""

if __name__ == "__main__":
    from numpy.lib.recfunctions import append_fields
    import sys
    if len(sys.argv) < 2:
        print(HELP_MESSAGE)
        exit()

    # Read in the timing data
    data = read_struct("IoZone-timetable-2017-10.csv")
    x = data[MODELLED_COLUMNS]
    y = data[TARGET_COLUMN]

    # Convert repeated tests into averages
    unique_x = set()
    all_y = {}
    for row,time in zip(x,y):
        row = tuple(row)
        unique_x.add( row )
        all_y[row] = all_y.get(row,[]) + [time]

    # Get the x and y in a form that can be modelled
    new_x = np.array([list(row) for row in unique_x], dtype=float)
    new_y = np.array([sum(all_y[row]) / len(all_y[row]) for row in unique_x])

    x_shift = np.min(new_x, axis=0)
    new_x -= x_shift
    x_scale = np.max(new_x, axis=0)
    new_x /= x_scale

    # Fit a delaunay model to the data
    model = Delaunay()
    model.fit(new_x, new_y)
    config_file = sys.argv[1]
    # Now read in the desired test settings
    test_points = read_struct(config_file)
    # Extract the columns that are provided
    test_cols = [c for c in MODELLED_COLUMNS+[ITERATION_COLUMN]
                 if c in test_points.dtype.names]
    test_points = test_points[test_cols].copy()
    # Fill in values for the columns that were not provided
    for c in [c for c in MODELLED_COLUMNS+[ITERATION_COLUMN]]:
        if c not in test_cols:
            print("Adding column of ones for '%s'."%(c))
            test_points = append_fields(
                test_points, c, np.ones(len(test_points)))
    # Extract out the number of iterations for each setting
    iterations = test_points[ITERATION_COLUMN].copy()
    test_points = test_points[MODELLED_COLUMNS].copy()
    if "DOE.Configuration.Table.csv" in config_file:
        # Set the number of iterations to 150
        iterations = 150
        # Scale the frequency into the same range as the training data
        test_points[MODELLED_COLUMNS[0]] *= 1000000
    # Get the settings into a standard numpy array
    test_points = np.array([list(r) for r in test_points[MODELLED_COLUMNS].copy()],dtype=float)
    # Normalize the new points according to the inputs
    test_points -= x_shift
    test_points /= x_scale
    # Guess times using the model
    guessed_times = model(test_points)
    # Multiply by iterations and sum to get the total estimate
    time_in_seconds = sum(guessed_times * iterations)
    # Convert total seconds estimate into days, hours, minutes.
    days = time_in_seconds // (60*60*24)
    time_in_seconds -= days * (60*60*24)
    hours = time_in_seconds // (60*60)
    time_in_seconds -= hours * (60*60)
    minutes = round(time_in_seconds / 60)
    output = []
    if days > 0:  output += [("%i day" + ("s" if days != 1 else ""))%(days)]
    if hours > 0: output += [("%i hour" + ("s" if hours != 1 else ""))%(hours)]
    if minutes > 0: output += [("%i minute" + ("s" if minutes != 1 else ""))%(minutes)]
    # Display output to user
    print(", ".join(output),"on one machine.")

    num_machines = 12
    # Multiply by iterations and sum to get the total estimate
    time_in_seconds = sum(guessed_times * iterations) / num_machines
    # Convert total seconds estimate into days, hours, minutes.
    days = time_in_seconds // (60*60*24)
    time_in_seconds -= days * (60*60*24)
    hours = time_in_seconds // (60*60)
    time_in_seconds -= hours * (60*60)
    minutes = round(time_in_seconds / 60)
    output = []
    if days > 0:  output += [("%i day" + ("s" if days != 1 else ""))%(days)]
    if hours > 0: output += [("%i hour" + ("s" if hours != 1 else ""))%(hours)]
    if minutes > 0: output += [("%i minute" + ("s" if minutes != 1 else ""))%(minutes)]
    # Display output to user
    print(", ".join(output),"on %i machines."%(num_machines))
    print(" at about %.1f seconds per test."%(sum(guessed_times) / len(test_points)))
