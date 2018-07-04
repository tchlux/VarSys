import os, random
from util.data import Data
from util.plot import Plot
from util.stats import cdf_fit_func
import numpy as np

# Should a subset of data be used
sub_data = True
# Source file and loaded data names.
_no_aesni = ""
# _no_aesni = "_no_aesni"
file_name = f"study{_no_aesni}.csv"
data_name = f"study{_no_aesni}.pkl.gz"
# sub_data_name = f"sub_study{_no_aesni}_10K.pkl"
sub_data_name = f"sub_study{_no_aesni}_1M.pkl"
sub_size = 10000 if "10K" in sub_data_name else 1000000
if sub_data:
    # Load data (and get a subset of it for testing and analysis)
    if not os.path.exists(sub_data_name):
        print("Loading file..")
        if not os.path.exists(data_name):
            d = Data.load(file_name, sep=",", verbose=True)
            print("Saving data..")
            d.save(data_name)
        else:
            d = Data.load(data_name)
        print("Done loading.")
        print(d)
        print()
        print("Getting random subsample..")
        sub_d = d[sorted(random.sample(range(len(d)), sub_size))]
        print("Saving random subsample..")
        sub_d.save(sub_data_name)
        d = sub_d
    else:
        print("Loading data from file..")
        d = Data.load(sub_data_name, verbose=True)
else:
    if not os.path.exists(data_name):
        print("Loading file..")
        d = Data.load(file_name, sep=",", verbose=True)
        print("Saving data..")
        d.save(data_name)
    else:
        print("Loading compressed data..")
        d = Data.load(data_name)
 
print("Getting the predictors..")
predictors = d.names[:2] + d.names[3:]
target = [d.names[2]]
print("Getting the unique values..")
unique_data = d[predictors].unique()
unique_data.sort()
print("Collecting by unique values..")
d = unique_data.collect(d)

print("Identifying rows that need to be removed..")
to_remove = [i for i in range(len(d)) if (len(d[i, d.names[-1]]) <= 1)]
print("Removing empty rows..")
for idx in to_remove[::-1]: d.pop(idx)

print("Generating distribution functions..")
d[d.names[-1] + " Distribution"] = map(cdf_fit_func, d[d.names[-1]])

print("Collecting system data into time models..")
system_settings = d.names[:2]
time_model_settings = d.names[2:4]
system_data = d[system_settings].unique()
time_models = []
for row in system_data:
    print(f"  {row}\r", end="", flush=True)
    time_data = d[d[system_settings[0]] == row[0]]
    time_data = time_data[time_data[system_settings[1]] == row[1]]
    time_models.append( time_data[time_data.names[2:]] )
system_data["Time Model"] = time_models
print("Saving processed data..")
system_data.save(f"processed_study{_no_aesni}.dill")

# Generate new data that is:
#   L3 | Packet | Time Model Data
# 
#  where "Time Model Data" is:
#   Byte | Value | Clock Cycles | CC CDF

