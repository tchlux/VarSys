from util.data import Struct
from util.stats import cdf_fit_func, ks_diff
import os



# Load random_readers data
if not os.path.exists("data_random_readers_[numeric].pkl.gz"):
    data = Struct().load("data_random_readers.pkl.gz")
    headers = [h for (h,t) in zip(data.names,data.types) 
               if (t in {int, float}) and (h not in {"Iteration","Runtime"})]
    num_data = data[headers]
    num_data.save("data_random_readers_[numeric].pkl.gz")
else:
    num_data = Struct().load("data_random_readers_[numeric].pkl.gz")

# Load distribution data
if not os.path.exists("data_random_readers_[cdfs].dill.gz"):
    # Reduce data into unique configurations.
    print(num_data)
    counts = {}
    for row in num_data:
        setting = tuple(row[:-1])
        counts[setting] = counts.get(setting, []) + [row[-1]]

    print("Unique configs: ", len(counts))
    print("Unique num runs:", sorted(set(len(counts[s]) for s in counts)))
    print("Less than 150:  ", sum(1 for s in counts if len(counts[s]) < 150))
    # Remove configurations with too little data
    for s in list(counts):
        if len(counts[s]) < 150: 
            print("Too little data:", s)
            counts.pop(s)
        else:
            counts[s] = counts[s][:150]

    # Process raw data into CDF functions.
    print()
    print("Generating CDFs for each configuration.")
    dist_struct = Struct(names=(num_data.names[:-1] + ["Throughput CDF"]))
    for s in sorted(counts):
        dist_struct.append( list(s) + [cdf_fit_func(counts[s])] )

    print()
    print(dist_struct)
    print()
    dist_struct.save("data_random_readers_[cdfs].dill.gz")
else:
    dist_struct = Struct().load("data_random_readers_[cdfs].dill.gz")


print()
print(dist_struct)
print()
data, info = dist_struct[dist_struct.names[:-1]].to_numpy_real()
print(data.shape)
print(data)
print()

from util.algorithms import Delaunay
import numpy as np

# num_batches = 4
# for batch in (0,):
ks_differences = []
for row in range(len(dist_struct)):
    # # Skip all but 1/3 of this data
    # if not ((len(dist_struct) / num_batches) * batch <= row 
    #         < (len(dist_struct) / num_batches) * (batch+1)):
    #     ks_differences.append(None)
    #     continue
    print(f"{row} : {len(dist_struct)}", end="\r")
    model = Delaunay()
    sub_data = np.vstack((data[:row], data[row+1:]))
    # Normalize the data to be in the unit hypercube
    shift = np.min(sub_data, axis=0)
    sub_data -= shift
    scale = np.max(sub_data, axis=0)
    sub_data /= scale
    # Fit a model to the data
    model.fit(sub_data)
    # Generate an estimated set of points and weights for this new spot
    pts, wts = model.points_and_weights((data[row] - shift) / scale)
    # Increment the points to align correctly with original indices
    pts = np.where(pts >= row, pts+1, pts)
    # Get the minimum and maximum value for all source CDFs
    min_max = np.array([dist_struct[p][-1]() for p in pts])
    min_max = (min(min_max[:,0]), max(min_max[:,1]))
    # Generate a guess CDF by taking a weighted sum of existing ones
    guess_cdf = lambda x=None: (
        sum(dist_struct[p][-1](x) * w  for (p,w) in zip(pts,wts))
        if (type(x) != type(None)) else min_max
    )
    # Calculate an error and record it
    error = ks_diff(dist_struct[row,-1], guess_cdf)[0]
    ks_differences.append(error)

# ks_differences = ks_differences + [None] * (len(dist_struct) - len(ks_differences))
print("Adding KS Stat error column")
dist_struct.add_column(ks_differences, name="KS Statistic for Prediction")
print("Constructing saveable struct")
to_save = dist_struct[dist_struct.names[:-2] + [dist_struct.names[-1]]]
print("Saving CSV")
to_save.save(f"prediction_results.csv")
