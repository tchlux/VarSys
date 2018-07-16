import os
import numpy as np
from util.data import Data

reduced_file = "time_model_mean_var.pkl"

if not os.path.exists(reduced_file):
    _no_aesni = ""
    file_name = f"processed_study{_no_aesni}.dill"

    print("Loading full file..")
    d = Data.load(file_name)
    print(d)
    # A function for overwriting the "Time Model" column
    def reduce_time_model(tm):
        tm.pop("Clock Cycles Distribution")
        tm["Mean Clock Cycles"] = map(lambda l: float(np.mean(l)), tm["Clock Cycles"])
        tm["Var Clock Cycles"] = map(lambda l: float(np.var(l)), tm["Clock Cycles"])
        tm.pop("Clock Cycles")
        # print(tm)
        # print(tm.to_matrix())
        # exit()
        return tm
    d["Time Model"] = map(reduce_time_model, d["Time Model"])
    print("Saving data..")
    d.save(reduced_file)
    print()
else:
    print("Loading data..")
    d = Data.load(reduced_file)
    print()



tm_names = d[0,-1].names
tm_types = d[0,-1].types
tm_shape = (len(d[0,-1]), len(d[0,-1].names) - 2)
tm_consts = d[0,-1].to_matrix().data[:,:2]
print(tm_names, tm_shape)
print(tm_consts)
truth = d.copy()
print()

# Convert data into numpy matrices
print("Converting data into prediction form..")
time_models = Data()
for row in d:
    time_models.append( list(map(float,row[-1].to_matrix().data[:,2:].flatten())) )
d.pop(d.names[-1])
d += time_models
print(len(d), len(d.names))
print(d.names[:10])
# d[d.names[-1]] = map(lambda tm: tm.to_matrix().data.flatten(), d[d.names[-1]])
print()

# Make predictions
from util.algorithms import Delaunay
print("Making predictions..")
predictions = d.predict(d.names[2:], model=Delaunay(), k=len(d))\
               .to_matrix().data.reshape((len(d), tm_shape[0], 2))
predictions = [Data(data=np.hstack((tm_consts, p)), names=tm_names, types=tm_types)
               for p in predictions]
truth["Predicted Time Models"] = predictions
print(truth)

tm_file = "time_model_predicted[Delaunay].pkl"
truth.save(tm_file)
