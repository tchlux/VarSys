import os
import numpy as np
from util.data import Data
from util.stats import cdf_fit
from util.algorithms import Delaunay

data_file_name = "[2018-10]_Xiaodong_Data.csv"
data_file = "data.pkl"

if not os.path.exists(data_file):
    d = Data.load(data_file_name, sep=",", verbose=True)
    print("Creating time model column..")
    d["Time Model"] = ((i,v,c) for i,v,c in zip(
        d["Byte index"], d["Byte possible value"], d["Number of Clock Cycles"]))
    print("Popping old columns..")
    d.pop("Byte index")
    d.pop("Byte possible value")
    d.pop("Number of Clock Cycles")
    print("Collecting time models by unique configuration..")
    d = d[d.names[:-1]].unique().collect(d)
    print("Organizing time models by index-value pairs..")
    d["Index-Values"] = (tuple(v[:-1] for v in tm) for tm in d["Time Model"])
    d["Time Model Size"] = map(len, d["Index-Values"])
    print("Converting time models into arrays..")
    d["Time Model"] = (np.array(sorted(tm)) for tm in d["Time Model"])
    print()
    print(d.shape)
    print(d.names)
    print()
    print("Saving data..")
    d.save("data.pkl")
else:
    d = Data.load(data_file)

# Get the bad configurations
bad_configs = Data()
unique_tms = {}
for i,tm in enumerate(d["Index-Values"]):
    unique_tms[tm] = unique_tms.get(tm, []) + [i]
for utm in sorted(unique_tms):
    if len(unique_tms[utm]) == 1:
        bad_configs += d[ unique_tms[utm][0], d.names[:5]+[-1] ]

print()
bad_indices = [unique_tms[utm][0] for utm in sorted(unique_tms) if len(unique_tms[utm])==1]
print("Bad configurations at indices", bad_indices)
print(bad_configs)

# d.summarize()
results_file = f"Time_Model_Predictions_Delaunay_{len(d)}-fold"

if not os.path.exists(results_file + ".pkl"):

    names = ["Private cache size", "Private cache associativity",
             "Shared cache size", "Shared cache associativity",
             "Packet size", "Byte index", "Byte possible value",
             "Predicted Number of Clock Cycles"]
    types = [int, int, int, int, int, int, int, float]
    results = Data(names=names, types=types)

    for i,(train, test) in enumerate(d.k_fold(k=len(d), only_indices=True)):
        print(f"Fold {i+1} of {len(d)}..", end="\r")
        # Remove the "bad" configurations from training data.
        for i in bad_indices:
            if i in train: train.remove(i)
        # Get the matrix form for the test data based on the training data.
        num = d[train, :5].to_matrix()
        # Construct the numeric training and testing data based on the training.
        train_x = num.data
        train_y = [tm for (i,tm) in enumerate(d["Time Model"]) if i in train]
        test_x = np.array([num.to_real(d[row][:5]) for row in test])
        test_y = [tm for (i,tm) in enumerate(d["Time Model"]) if i in test]
        # Normalize the training and testing coordinates.
        train_x = (train_x - num.shift) / num.scale
        test_x = (train_x - num.shift) / num.scale
        # Make a prediction with model.
        model = Delaunay()
        model.fit(train_x, train_y)
        guess_y = model(test_x)
        # Store the results in a data object.
        config = d[test[0]][:5]
        for (index, value, clock) in guess_y[0]:
            index = int(round(index))
            value = int(round(value))
            results.append( config + [index, value, clock] )

    results.save(results_file+".pkl")
    results.save(results_file+".csv.gz")

else:
    results = Data.load(results_file+".pkl")

print("data.shape:    ",'(272996, 8)')
print("results.shape: ",results.shape)


# Collect data by all the settings that would be changed in reality
# (to affect the "Time Model") 
# -> Private cache size
#    Private cache associativity
#    Shared cache size
#    Shared associativity
#    Packet size

# The "Time model" is a model showing how long different byt positions
# and values take to execute (in clock cycles):
# -> Number of clock cycles
#    Byte index
#    Byte possible value

# The latter two parameters form a grid of values, at each location
# there is a distribution of numbers of clock cycles.
