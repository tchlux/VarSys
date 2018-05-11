import numpy as np
from util.data import read_struct
import pickle

data = read_struct("stage1.csv")
print(data.dtype.names)
print(data.shape)
# Reduce to interesting columns
cols_to_keep = ["Frequency", "Record_Size", "Threads", "Throughput", "Test_Type"]
data = data[cols_to_keep]
# Only look at "Fread"
data = data[data["Test_Type"] == "freaders"]
# Remove the "Test_Type" column entirely
data = data[cols_to_keep[:-1]]
# Save the data pickle to file
with open("Stage_1_Fread.pkl", "wb") as f:
    pickle.dump(data, f)
