import os
import numpy as np
from util.data import Data
# from util.algorithms import Delaunay as Interpolant
# from util.algorithms import Voronoi as Interpolant
from util.algorithms import NearestNeighbor as Interpolant
# from sklearn.neural_network import MLPRegressor as Regressor
from sklearn.metrics import confusion_matrix

# =========================
#      Process Li Data     
# =========================

if not os.path.exists("dist_predicted.pkl"):
    print("Reading data from file..")
    raw_data_file = "EM_results.csv"
    d = Data.load(raw_data_file)
    config = d.names[1:6]
    predict = d.names[6:]
    d += d.predict(predict, model=Regressor())
    d.save("dist_predicted.pkl")
else:
    d = Data.load("dist_predicted.pkl")

outcomes = d[d.names[-4:]]
# Measure performance of type-of-distribution prediction

labels = sorted(set(outcomes['disname']))
perf = confusion_matrix(list(outcomes["disname"]), 
                        list(outcomes["disname Predicted"]),
                        labels=labels)
print()
print(perf.shape)
print("Classification Contigency Table:")
print(labels)
print(perf)
print()
print("Correct classifications:\n", np.sum(np.diag(perf)) / np.sum(perf))
print()

# Measure performance of number-of-distributions prediction
errors = (np.array(list(outcomes["disnumcomp Predicted"])) -
          np.array(list(outcomes["disnumcomp"])))
correct = sum(errors == 0) / len(errors)
off_by_1 = sum(abs(errors) == 1) / len(errors)
off_by_2 = sum(abs(errors) == 2) / len(errors)
print(f"Correct disnumcomp predictions: {correct:.2f}")
print(f"Absolute error == 1:            {off_by_1:.2f} ({correct + off_by_1:.2f})")
print(f"Absolute error == 2:            {off_by_2:.2f} ({correct+off_by_1+off_by_2:.2f})")

d.save("EM_Prediction.csv")

exit()

from util.plotly import Plot
p = Plot("Guess vs. Truth (number of distributions)", "Guess", "Truth")
p.add("Data",list(outcomes["predicted disnumcomp"]), 
      list(outcomes["disnumcomp"]), show_in_legend=False,
      color=p.color(1))
p.show()
# 
p = Plot("Histogra of errors when prediction number of distributions")
p.add_histogram("Data", errors, color=p.color(1), show_in_legend=False)
p.show(append=True)
exit()

# - send Li my prediction results
# - send Amy email about video
