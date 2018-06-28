import os
from util.plot import Plot
from util.data import Data
from util.system import pause

# File names for demo purposes
home_dir = os.path.expanduser("~")
data_file_name = os.path.join(home_dir,"Git","VarSys",
                              "Data_[2017-01-13]_THOMAS_RESULTS.csv")

# Load the data from file
print("Loading data..")
data = Data.load(data_file_name, sep=",")
print(data)
print()
pause()

# Print out a summary of the data
data.summarize()
pause()

# Predict the values for one of the columns
print("Running predictions..")
data += data.predict(data.names[-1])

# Get column names
print()
truth_col = data[data.names[-2]]
guess_col = data[data.names[-1]]
error_col_name = data.names[-1] + " Error"
print("truth col name: ",data.names[-2])
print("guess col name: ",data.names[-1])
print("error col name: ",error_col_name)
print()
pause()

# Add a column for the error (using a generator) and print updated data
print("Computing errors..")
data[error_col_name] = ((g-t) for (t,g) in zip(truth_col, guess_col))
print(data)
pause()

# Generate a histogram of the error.
print("Generating histogram of errors..")
p = Plot()
p.add_histogram(error_col_name, data[error_col_name])
p.plot(file_name="demo_prediction_errors.html", show=True, 
       show_legend=False)
pause()

# Get a slice of the data (for 3D plotting)
print("Reducing and organizing plot data..")
X = "F_size"
Y = "R_Size"
Z = "Throughput"
# Reduce the data to only the "read" mode
data = data[ data["Mode"] == "Read" ]
plot_data = data[X,Y,Z,Z+" Predicted"]
pause()

# Generate a 3D plot of the truth values and the predictions.
p = Plot()
p.add("Known Throughputs", plot_data["F_size"], plot_data["R_Size"],
      plot_data["Throughput"])
p.add("Predicted Throughputs", plot_data["F_size"], plot_data["R_Size"],
      plot_data["Throughput Predicted"])
p.plot(file_name="demo_throughput_3D.html", show=True)
