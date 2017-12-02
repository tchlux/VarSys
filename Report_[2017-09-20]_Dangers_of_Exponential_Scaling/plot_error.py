from util.plotly import Plot
from util.data import read_struct
import numpy as np


data = read_struct("New_Predictions.csv")

# Extract those points for testing and the repeated values
testing_f_r_sizes = [(512,32),(512,128),(512,256),(768,32),(768,128)]
numeric_columns = ["File_Size","Record_Size"]

# Reduce testing data to those specific rows
testing_indices = [tuple(row) in testing_f_r_sizes
                   for row in data[numeric_columns]]
# Extract the settings for known values from the test data
repeat_indices = [not v for v in testing_indices]
repeat_data = data[repeat_indices]
data = data[testing_indices]

# Headers and some targets
header = list(data.dtype.names)
truth = "True_Std"
original = "Original_Std"
alg_headers = header[header.index(truth)+1:]

p = Plot("Distribution of Errors Per Algorithm")

low = -100 * 10**6
high = 75 * 10**6
size = (high - low) / 100
bins = dict(start=low, end=high, size=size)
print(bins)

errors = {}

for a in alg_headers:
    errors[a] = (data[truth], data[a])

true_error = (repeat_data[truth], repeat_data[original])
errors[original] = true_error

# Plot the errors

def error_1(truth, pred):
    total = 0.0
    for t,p in zip(truth, pred):
        total += abs(p/t - 1)
    return total / len(truth)

def error_2(truth, pred):
    numerator = 0.0
    denominator = 0.0
    for t,p in zip(truth, pred):
        numerator += (p - t)**2
        denominator += t
    numerator = (numerator / len(truth)) ** (0.5)
    denominator /= len(truth)
    return numerator / denominator

def error_3(truth, pred):
    numerator = 0.0
    denominator = 0.0
    for t,p in zip(truth, pred):
        numerator += abs(p - t)
        denominator += t
    return numerator / denominator

names = ["Error 1", "Error 2", "Error 3"]
statistics = [error_1, error_2, error_3]


for n,error_measure in zip(names,statistics):
    results = []
    for algorithm in sorted(errors.keys()):
        truth, pred = errors[algorithm]
        error = error_measure(truth, pred)
        results.append( (algorithm, error) )
    results.sort( key=lambda i: i[1])
    print(n)
    for v in results:
        print(v)
    print()
    print()

exit()


#      Plotting the historgrams     
# ==================================
for a in sorted(errors.keys()):
    p.add_histogram( a, errors[a][0] - errors[a][1], opacity=0.6
                     ,autobinx=False, xbins=bins )
p.plot(x_range=[low,high])
