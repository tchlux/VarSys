from util.data import Data

print("Loading original data..")
d = Data.load("2019-03-19_Xiaodong-data.csv")
# Collect by unique configurations.
print("Collecting data by four parameters..")
d = d[:,:-1].unique().collect(d)
d["Samples"] = map(len, d["Clock Cycles"])
d["Mean Clock Cycles"] = map(lambda l: sum(l) / len(l), d["Clock Cycles"])
d.reorder(d.names[:-3] + ["Mean Clock Cycles", "Samples"])
d["Index"] = range(len(d))
d["Predicted Mean Clock Cycles"] = (0. for row in d)
# Remove the unused columns.
d.pop("Samples")
d.pop("Clock Cycles")
# Print the data for verification.
print(d)
print()
print("Making predictions..")


from util.approximate import Delaunay as Algorithm
import numpy as np
for i,(train, test) in enumerate(d.k_fold(k=len(d))):
    print(f" fold {i+1}..", train.shape, test.shape)
    matrix = train[:,:4].to_matrix()
    train_x = matrix.data
    train_y = np.array(list(train["Mean Clock Cycles"]))
    # Assume leave-one-out prediction.
    test_x = np.array(test[0][:4], dtype=float)
    # normalize the train and test by the train
    train_x -= matrix.shift
    train_x /= matrix.scale
    test_x -= matrix.shift
    test_x /= matrix.scale
    model = Algorithm()
    model.fit(train_x, train_y)
    # Store the output predictions.
    d[test[0,"Index"], "Predicted Mean Clock Cycles"] = float(model.predict(test_x))

# Remove the additional index column that was added.
d.pop("Index")
# Compute the prediction errors.
d["Error"] = (g-t for (t,g) in zip(d["Mean Clock Cycles"], d["Predicted Mean Clock Cycles"]))
d["Absolute Error"] = (abs(g-t) for (t,g) in zip(d["Mean Clock Cycles"], d["Predicted Mean Clock Cycles"]))
print()
print(d)

d.summarize()

d.save("predictions.csv")
