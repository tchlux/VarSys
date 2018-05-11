import os
from util.data import Struct
from util.multi_dim_analysis import get_rand_sets_of_points
from util.plotly import Plot
from setup import MaxBoxMesh, IterativeBoxMesh, VoronoiMesh
from util.algorithms import Delaunay
import numpy as np

# ===============================
#      Load Performance Data     
# ===============================

performance_data_file = "performance_data.pkl"
if not os.path.exists(performance_data_file):
    # For processing raw data
    readers_data_file = "readers_data.pkl"
    if not os.path.exists(readers_data_file):
        print("Loading data from file..")
        data = Struct.load("data_full.pkl.gz")
        interest_cols = ['Frequency', 'File Size', 'Record Size', 
                         'Num Threads', 'Throughput']

        print("Reducing data to 'readers' test..")
        reduced_data = data[data["Test"] == "readers"]
        del data

        print("Reducing data to columns of interest..")
        data = reduced_data[interest_cols]
        del reduced_data

        print(data)
        data.save(readers_data_file)
    else:
        print(f"Loading {readers_data_file}..")
        data = Struct.load(readers_data_file)

    # Reduce the data to only those unique configuations (and their performances)
    print("Reducing to unique configurations..")
    unique_configs = {config : [] for config in map(tuple,data[data.names[:-1]])}
    for row in data:
        config = tuple(row[:-1])
        unique_configs[config] += [row[-1]]
    new_names = data.names[:-1] + ["Throughputs"]
    del data

    # Clean up the performance data and save to file.
    print("Cleaning up data and converting to Struct..")
    performance_data = Struct(names=new_names)
    for config in unique_configs:
        if (len(unique_configs[config]) < 150):
            print(f"{len(unique_configs[config])} samples for setting {config}, throwing out configuration..")
            continue
        performance_data.append( config + (unique_configs[config][:150],) )
    del unique_configs

    # 145 samples for setting (2300000, 1024, 128, 16), throwing out configuration..
    # 100 samples for setting (2300000, 1024, 128, 24), throwing out configuration..
    # 59 samples for setting (2300000, 1024, 128, 32), throwing out configuration..
    # 31 samples for setting (2300000, 1024, 128, 40), throwing out configuration..

    print("Adding variance column..")
    performance_data.add_column(
        map(lambda p: float(np.var(p)), performance_data["Throughputs"]),
        "Variance"
    )

    print("Adding mean column..")
    performance_data.add_column(
        map(lambda p: float(np.mean(p)), performance_data["Throughputs"]),
        "Mean"
    )

    print("Saving performance data..")
    performance_data.save(performance_data_file)
else:
    print("Loading performance data..")
    performance_data = Struct.load(performance_data_file)


all_data = performance_data[performance_data.names[:-3] + ["Variance"]]

# Get the data in a matrix
data, info = all_data.to_numpy_real()

# Normalize the data independent values
shift = np.min(data[:,:-1], axis=0)
data[:,:-1] -= shift
scale = np.max(data[:,:-1], axis=0)
data[:,:-1] /= scale

# Suffle the data
np.random.seed(0)
np.random.shuffle(data)
# Reduce data to 1000 points
data = data[:1000]

# Get the number of training points
num_train = int(.5+.8*len(data))

# Get the training points
train_points_file = "train_points.py"
if not os.path.exists(train_points_file):
    help(get_rand_sets_of_points)
    print("Creating well-spaced training points..")
    train_points = get_rand_sets_of_points(list(range(len(data))), num_train,
                                    num_sets=1, data=data[:,:-1],
                                    well_spaced=True)[0]
    with open(train_points_file, "w") as f:
        print(f"train_points = {train_points}", file=f)
else:
    print("Loading training points..")
    from train_points import train_points
# Get the testing points
test_points = [i for i in range(len(data)) if i not in train_points]


train_x, train_y = data[train_points,:-1], data[train_points,-1]
test_x, test_y = data[test_points,:-1], data[test_points,-1]

print(train_x.shape)
print(test_x.shape)

algorithms = [IterativeBoxMesh, VoronoiMesh, Delaunay]
error_tols = [.67, .4, None]
for alg,et in zip(algorithms, error_tols):
    # Handle the two different types of initialization
    if (alg == Delaunay):
        model = alg()
        model.control_points = np.ones((1,1))
        model.coefficients = model.control_points
    else:
        model = alg(et)
    print()
    print(f"Fitting {alg.__name__}..")
    # Get the LS fit guesses
    model.fit(train_x.copy(), train_y.copy())
    print(f"Using {len(model.control_points[0])} control points..")
    print(f"Evaluating {alg.__name__} for LS..")
    ls_guesses = model(test_x.copy())
    
    # Update coefficients to make this a lagrange interpolant
    print("Updating coefficients..")
    for i in range(len(model.control_points)):
        closest_idx = np.argmin(np.sum((train_x - model.control_points[:,i])**2,axis=1))
        model.coefficients[i] = train_y[closest_idx]
    print(f"Evaluating {alg.__name__} for lagrangian..")
    lg_guesses = model(test_x)

    # print("Plotting results..")
    # # Plot the errors for the algorithm using LS versus Lagrange interpolation
    # p = Plot(f"Signed Relative Errors for {alg.__name__}", "Error", "Count")
    # p.add_histogram(f"Least Squares {alg.__name__}", (test_y - ls_guesses)/test_y,
    #                 start=-500, end=100)
    # p.add_histogram(f"Lagrange      {alg.__name__}", (test_y - lg_guesses)/test_y,
    #                 start=-500, end=100)
    # p.plot(x_range=[-500,100], file_name="ls_vs_lg.html", append=(alg!=algorithms[0]))
    print("                     LS Fit       Lagrange")
    for perc in [0, 25, 50, 75, 100]:
        print(f"{perc: 5d}th percentile: {np.percentile(abs(ls_guesses),perc): 5.4e}  {np.percentile(abs(lg_guesses),perc): 5.4e}")
