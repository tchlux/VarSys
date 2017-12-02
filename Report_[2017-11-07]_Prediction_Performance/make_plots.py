from util.data import read_struct
from util.plotly import *
import pickle
from util.data import read_struct

# DATA_FILE = "Iozone-2016.pkl"
# DATA_FILE = "Data.pkl"

# with open(DATA_FILE, "rb") as f:
#     data = pickle.load(f)

data = read_struct("HPC_Compare_Set.csv")

#  ALGORITHM     MEAN ABSOLUTE ERROR

# 1 Dimension                             # 1 Dimension                        
#  Delaunay      395400                   #  Delaunay        2060267           
#  MaxBoxMesh    399686                   #  FinalDelaunay   2060267           
#  LSHEP         516740                   #  NewDelaunay     2060267           
#  MARS         1050553                   #  OldDelaunay     2060267           
                                          #  qHullDelaunay   2060267           

# 2 Dimensions                            # 2 Dimensions                                                            
#  LSHEP         283206                   #  qHullDelaunay   2877587           
#  Delaunay      405099                   #  NewDelaunay     2878570           
#  MaxBoxMesh    453502                   #  Delaunay        2883286           
#  MARS          605142                   #  FinalDelaunay   2914714           
                                          #  OldDelaunay     3566723           
                                          
# ALGORITHM   AVG ABSOLUTE MEAN ERROR     # ALGORITHM   AVG ABSOLUTE MEAN ERROR
# LSHEP         399973                    # qHullDelaunay   2468927            
# Delaunay      400250                    # NewDelaunay     2469418            
# MaxBoxMesh    426594                    # Delaunay        2471777            
# MARS          827848                    # FinalDelaunay   2487490            
                                          # OldDelaunay     2813495            

print(data.dtype.names)
print(len(data))

algorithms = sorted(np.unique(data["Algorithm"]))

# lower = -5000000
# upper =  5000000
lower = -2
upper =  2
y_bins = dict(start=lower, end=upper, size=(upper - lower) / 100)
lower_x = 0
upper_x = .21

# ============================================
#      Making histogram plot by dimension     
# ============================================

results = {}

plots = []
out_name = "Prediction-Performance.html"
# out_name = "Prediction-Performance-Histograms-by-Dimension.html"
brief_title = "Distributions of Prediction Error with Increasing Dimension Data"
full_text_title = "<b>"+brief_title+"</b><br>"+\
                  "Where the individual x-axes are the probability ranges"+\
                  " for the bars in each histogram."
error_col = "Relative_Mean_Error"

for dim in range(1,np.max(data["Dimension"])+1):
    dim_data = data[(data["Dimension"] == dim)]
    title = (brief_title if dim == 4 else "")
    x_axis = "%i Dimension"%dim + ("s" if dim > 1 else "")
    y_axis = ("Relative Error in Predicted System Throughput" if dim == 1 else "")
    p = Plot(title, x_axis, y_axis)
    for alg in algorithms:
        alg_data = dim_data[(dim_data["Algorithm"] == alg)]
        values = alg_data[error_col]
        y_bins['start'] = min(values)
        y_bins['end'] = max(values)
        results[((dim,),alg)] = abs(values).mean()
        p.add_histogram(alg, values, bar_spacing="y", group=alg,
                        histnorm='probability',
                        ybins=y_bins, show_in_legend=(dim==1))
    plots.append(p.plot(y_range=(lower,upper), x_range=(lower_x,upper_x),
                         html=False, layout={'title':full_text_title}))

# multiplot(plots, shared_y=True, file_name=out_name,show=False)

# ======================================================
#      Making histogram plot by training percentage     
# ======================================================

plots = []
# out_name = "Prediction-Performance-Histograms-by-Train-Percentage.html"
brief_title = "Distributions of Prediction Error with Increasing Amounts of Training Data"
full_text_title = "<b>"+brief_title+"</b><br>"+\
                  "Where the remaining data is used for testing the"+\
                  " predictive performance of the algorithms.<br>" +\
                  "The individual x-axes are the probability ranges"+\
                  " for the bars in each histogram."
training_ranges = [(0,25),(25,50),(50,75),(75,100)]
for min_train,max_train in training_ranges:
    # Some booleans used for naming
    first_entry = (min_train,max_train) == training_ranges[0]
    last_entry = (min_train,max_train) == training_ranges[-1]
    # Reduce to a local set of data
    set_data = data[(min_train <= data["Train_Percentage"])]
    set_data = set_data[(set_data["Train_Percentage"] < max_train)]
    # Generate the title and axis labels
    title = (brief_title if last_entry else "")
    min_val = min(set_data["Train_Percentage"])
    max_val = max(set_data["Train_Percentage"])
    x_axis = "[%.0f - %.0f%%] Training"%(min_val,max_val)
    y_axis = ("Relative Error in Predicted System Throughput" if first_entry else "")
    # Initialize thep lot
    p = Plot(title, x_axis, y_axis)
    for alg in algorithms:
        alg_data = set_data[(set_data["Algorithm"] == alg)]
        values = alg_data[error_col]
        y_bins['start'] = min(values)
        y_bins['end'] = max(values)
        results[((round(min_val),round(max_val)),alg)] = abs(values).mean()
        p.add_histogram(alg, values, bar_spacing="y", group=alg,
                        histnorm='probability', ybins=y_bins,
                        show_in_legend=first_entry)
    plots.append(p.plot(y_range=(lower,upper), x_range=(lower_x,upper_x),
                         html=False, layout={'title':full_text_title}))

# multiplot(plots, shared_y=True, file_name=out_name, append=True, show=True)

# ===========================================
#      Printing out results in text form     
# ===========================================

longest_name = max(map(len,algorithms))
longest_result = max(map(lambda v: len(str(round(v))), results.values()))
keys = sorted(results.keys(), key=lambda v: (v[0], results[v]))
print()
print(" ALGORITHM     MEAN ABSOLUTE ERROR")
first = None
alg_results = {a:[] for a in algorithms}
for k in keys:
    # Print headers for the various sections
    if k[0] != first:
        print()
        if len(k[0]) == 1:
            dim = k[0][0]
            print(dim, "Dimension" + ("s" if dim > 1 else ""))
        else:
            min_val, max_val = k[0]
            print("[%.0f - %.0f%%] Training Data"%(min_val,max_val))
    # Extract out information for printing neatly
    alg = k[1]
    val = "%0.2f"%(results[k])
    print('',alg+(longest_name - len(alg))*" ",
          (longest_result - len(val))*" " + val)
    # Collect algorithm results
    alg_results[alg].append(results[k])
    # Reset new loop check "first"
    first = k[0]
print()
print()
algorithms.sort(key=lambda alg: sum(alg_results[alg]) / len(alg_results[alg]))
print("ALGORITHM   AVG RELATIVE MEAN ERROR")
for alg in algorithms:
    val = "%0.2f"%(sum(alg_results[alg]) / len(alg_results[alg]))
    print(alg+(longest_name-len(alg))*" ",
          " "*(longest_result - len(val))+val)
exit()
