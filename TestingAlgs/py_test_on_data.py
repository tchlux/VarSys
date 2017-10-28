from util.plotly import Plot, multiplot
from py_settings import *
import numpy as np
import time

# MANUALLY TESTING FAILURE CASE
# train = "/Users/thomaslux/Desktop/TestingAlgs/VarSys-2016-Mean_MDA_Data_0/3-F_Size--R_Size--Freq/10.42--89.58/Train-0.csv"
# test = "/Users/thomaslux/Desktop/TestingAlgs/VarSys-2016-Mean_MDA_Data_0/3-F_Size--R_Size--Freq/10.42--89.58/Test-0.csv"
# # ^^ FAILURE CASE FOR MARS

# i = 150
# train = "/Users/thomaslux/Desktop/TestingAlgs/VarSys-2016-Mean_MDA_Data_0/2-Threads--Freq/50--50/Train-%s.csv"%i
# test = "/Users/thomaslux/Desktop/TestingAlgs/VarSys-2016-Mean_MDA_Data_0/2-Threads--Freq/50--50/Test-%s.csv"%i

# train = "/Users/thomaslux/Desktop/TestingAlgs/VarSys-2016-Mean_MDA_Data_0/4-F_Size--R_Size--Threads--Freq/25--75/Train-169.csv"
# test = "/Users/thomaslux/Desktop/TestingAlgs/VarSys-2016-Mean_MDA_Data_0/4-F_Size--R_Size--Threads--Freq/25--75/Test-169.csv"

# train = "/Users/thomaslux/Desktop/TestingAlgs/VarSys-2016-Mean_MDA_Data_0/4-F_Size--R_Size--Threads--Freq/20.02--79.98/Train-7.csv"
# test = "/Users/thomaslux/Desktop/TestingAlgs/VarSys-2016-Mean_MDA_Data_0/4-F_Size--R_Size--Threads--Freq/20.02--79.98/Test-7.csv"

# train = "/Users/thomaslux/Desktop/TestingAlgs/VarSys-2016-Mean_MDA_Data_0/4-F_Size--R_Size--Threads--Freq/9.95--90.05/Train-134.csv"
# test = "/Users/thomaslux/Desktop/TestingAlgs/VarSys-2016-Mean_MDA_Data_0/4-F_Size--R_Size--Threads--Freq/9.95--90.05/Test-134.csv"

# i = 3
# train = "/Users/thomaslux/Desktop/TestingAlgs/VarSys-2016-Mean_MDA_Data_0/3-R_Size--Threads--Freq/5.09--94.91/Train-1.csv"
# test = "/Users/thomaslux/Desktop/TestingAlgs/VarSys-2016-Mean_MDA_Data_0/3-R_Size--Threads--Freq/5.09--94.91/Test-1.csv"

# i = 51
# train = "/Users/thomaslux/Desktop/TestingAlgs/VarSys-2016-Mean_MDA_Data_0/2-Threads--Freq/4.86--95.14/Train-%i.csv"%i
# test = "/Users/thomaslux/Desktop/TestingAlgs/VarSys-2016-Mean_MDA_Data_0/2-Threads--Freq/4.86--95.14/Test-%i.csv"%i

# train = "/Users/thomaslux/Desktop/VarSys-2016-Mean_MDA_Data_0/2-Threads--Freq/4.86--95.14/Train-178.csv"
# test =  "/Users/thomaslux/Desktop/VarSys-2016-Mean_MDA_Data_0/2-Threads--Freq/4.86--95.14/Test-178.csv"

# train = "/Users/thomaslux/Desktop/VarSys-2016-Mean_MDA_Data_0/4-F_Size--R_Size--Threads--Freq/29.98--70.02/Train-135.csv"
# test =  "/Users/thomaslux/Desktop/VarSys-2016-Mean_MDA_Data_0/4-F_Size--R_Size--Threads--Freq/29.98--70.02/Test-135.csv"

train = "/Users/thomaslux/Desktop/VarSys-2016-Mean_MDA_Data_0/4-F_Size--R_Size--Threads--Freq/4.98--95.02/Train-175.csv"
test  = "/Users/thomaslux/Desktop/VarSys-2016-Mean_MDA_Data_0/4-F_Size--R_Size--Threads--Freq/4.98--95.02/Test-175.csv"

   # /Users/thomaslux/Desktop/VarSys-2016-Mean_MDA_Data_0/4-F_Size--R_Size--Threads--Freq/4.98--95.02/Train-175.csv
   # /Users/thomaslux/Desktop/VarSys-2016-Mean_MDA_Data_0/4-F_Size--R_Size--Threads--Freq/4.98--95.02/Test-175.csv

# train = "/Users/thomaslux/Desktop/VarSys-2016-Mean_MDA_Data_0/2-F_Size--Threads/14.81--85.19/Train-0.csv"
# test  = "/Users/thomaslux/Desktop/VarSys-2016-Mean_MDA_Data_0/2-F_Size--Threads/14.81--85.19/Test-0.csv"


# Create surface
# surf = NearestNeighbor()
# surf = MARS()
# surf = LSHEP()
surf = Delaunay()
# surf = MaxBoxMesh()
# surf = qHullDelaunay()


train_data = np.loadtxt(train, dtype=float, delimiter=",", ndmin=2)
test_data = np.loadtxt(test, dtype=float, delimiter=",", ndmin=2)
# Extract the points and values from the data file
train_points = train_data[:,:-1]
train_values = train_data[:,-1]
test_points = test_data[:,:-1]
test_values = test_data[:,-1]
# Normalize the data to the training-unit-hypercube
shift = train_points.min(axis=0)
train_points -= shift
test_points -= shift
scale = train_points.max(axis=0)
train_points /= scale
test_points /= scale



print()
print("Training data")
print("  Number    %i"%(train_values.shape[0]))
print("  Dim       %i"%(train_points.shape[1]))
print("  Pt Min    %s"%(train_points.min(axis=0)))
print("  Pt Max    %s"%(train_points.max(axis=0)))
print("  Pt Mean   %s"%(train_points.mean(axis=0)))
print("  Val Min   %.2f"%(train_values.min(axis=0)))
print("  Val Max   %.2f"%(train_values.max(axis=0)))
print("  Val Mean  %.2f"%(train_values.mean(axis=0)))
print()
print("Testing data")
print("  Number    %i"%(test_values.shape[0]))
print("  Dim       %i"%(test_points.shape[1]))
print("  Pt Min    %s"%(test_points.min(axis=0)))
print("  Pt Max    %s"%(test_points.max(axis=0)))
print("  Pt Mean   %s"%(test_points.mean(axis=0)))
print("  Val Min   %.2f"%(test_values.min(axis=0)))
print("  Val Max   %.2f"%(test_values.max(axis=0)))
print("  Val Mean  %.2f"%(test_values.mean(axis=0)))
print()
print()
# Make sure the entire array prints
np.set_printoptions(threshold=float('inf'))
clean_array = lambda arr: " "+str(arr).replace(",","").replace("[","").replace("]","")
print(clean_array(train_points))
print()

name = (str(surf.__class__).split("'")[1]).split(".")[1]
print()
print("Fitting...")
surf.fit(train_points, train_values)

print("Evaluating at test points...")
start = time.time()
test_predictions = surf(test_points)
print("  %.2f seconds"%(time.time()-start))

print("Evaluating at train points...")
start = time.time()
train_predictions = surf(train_points)
print("  %.2f seconds"%(time.time()-start))

# ========================================
#      CUSTOM CODE FOR DELAUNAY CHECK     
# ========================================

surf2 = qHullDelaunay()
surf2.fit(train_points, train_values)
test_predictions2 = surf2(test_points)
bad_ind = np.argmax(abs(test_predictions2 - test_predictions))

print()
print("="*70)
print("="*70)
print("="*70)

# bad_ind = 100
print("Bad index:", bad_ind)

for ind in [bad_ind]:
    # for ind in range(0,len(test_points)):
    print()
    print()
    print("Test point:")
    print(clean_array(test_points[bad_ind]))
    print()
    surf(test_points[bad_ind], debug=True, verbose=True)
    print()
    surf2(test_points[bad_ind], debug=True, verbose=True)
    print()

from util.notifications import send_email
del_simp = surf(test_points[bad_ind], debug=True)
qhull_simp = surf2(test_points[bad_ind], debug=True)
subject = "Delaunay Case of Interest"
body = "\nTrain Points:\n" + clean_array(train_points)
body += "\n\nTest Point:\n" + clean_array(test_points[bad_ind])
body += "\n\nDelaunay Simplex:  (0-indexed)\n" + clean_array(del_simp)
body += "\n\nQHull Simplex:     (0-indexed)\n" + clean_array(qhull_simp)
body += "\n\n~Thomas"

# print()
# print(body, end="")

# send_email(subject=subject, body=body, recipient="thchang@vt.edu", verbose=True)
exit()

# ========================================
# ========================================

test_error = abs(test_values - test_predictions) / test_values
print()
print("Error at test points:")
print("  Min    %.2f"%(test_error.min()))
print("  Max    %.2f"%(test_error.max()))
print("  Mean   %.2f"%(test_error.mean()))
print("  Stdev  %.2f"%(test_error.std()))

train_error = abs(train_values - train_predictions) / train_values
print()
print("Error at train points:")
print("  Min    %.2f"%(train_error.min()))
print("  Max    %.2f"%(train_error.max()))
print("  Mean   %.2f"%(train_error.mean()))
print("  Stdev  %.2f"%(train_error.std()))
print()

if test_points.shape[0] == 1:
    test_predictions = [test_predictions]

# Creating 3D visual of the surface over the data
if test_points.shape[1] <= 2:
    p = Plot()
    p.add("Train Points",*(train_points.T),train_values)
    p.add("Test Points",*(test_points.T),test_values)
    p.add(name,*(train_points.T),train_values,plot_type='surface', use_gradient=True)
    p.plot()

elif test_points.shape[1] == 3:
    p = Plot()
    p.add("Train Points",*(train_points.T))
    p.add("Test Points",*(test_points.T))
    p.plot()

