import os
from util.data import read_struct
from util.paper import latex_table


predictors = ["Mean", "Variance"]
algorithms = ["LSHEP", "MARS", "MLPRegressor", "qHullDelaunay", "SVR"]

# Process all raw files into pkl files that hold the contents (for
# easier loading later on)
for p in predictors:
    summary_output = p+"_all_summary.pkl"
    all_data = None
    for a in algorithms:
        source = "Output-%s_%s-MDA_results.csv"%(a,p)
        a_data = read_struct(source)
        print(len(a_data), source)
        
# Read in the data files, merge them into one pkl
# Process pkl as I do normally 
