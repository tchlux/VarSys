import os
import numpy as np
from util.data import read_struct
from util.plotly import Plot, multiplot
from util.algorithms import Delaunay, qHullDelaunay
from util import CLEAN_ARRAY_STRING


# Function for comparing the performance of the different delaunay
# techniques on a given training and testing file
def test_performance(train, test, plot=False):
    s1 = Delaunay()
    s2 = qHullDelaunay()
    # Get (and normalize)the training points
    train_points = np.loadtxt(train,delimiter=",",ndmin=2)
    train_values = train_points[:,-1]
    train_points = train_points[:,:-1]
    shift = np.min(train_points, axis=0)
    train_points -= shift
    scale = np.max(train_points, axis=0)
    train_points /= scale
    # Get (and normalize) the testing points
    test_points = np.loadtxt(test,delimiter=",",ndmin=2)
    test_values = test_points[:,-1]
    test_points = test_points[:,:-1]
    # normalize the test points
    test_points -= shift
    test_points /= scale
    # Fit the data
    s1.fit(train_points, train_values)
    s2.fit(train_points, train_values)
    # Calculate predicitons
    pred_1 = s1(test_points)
    pred_2 = s2(test_points)
    # Compute the differences between the two algorithms
    diff = abs(pred_1 - pred_2)
    bad_ind = np.argmax(diff)
    avg = (pred_1[bad_ind] + pred_2[bad_ind]) / 2
    perc_diff = diff[bad_ind] / avg
    # Plot if the user wants
    if plot:
        # Get the simplices returned by each delaunay algorithm
        simp1 = np.asarray(s1(test_points[bad_ind], debug=True), dtype=int)
        simp2 = np.asarray(s2(test_points[bad_ind], debug=True), dtype=int)
        if (all(s in simp1 for s in simp2) and all(s in simp2 for s in simp1)):
            return
        # Print the file names
        print()
        # Print the training points
        print(CLEAN_ARRAY_STRING(train_points))
        # Print the percent difference
        print()
        print("","Bad point:    ", (test_points[bad_ind] + shift) * scale)
        print("","Delaunay pred:", pred_1[bad_ind])
        print("","qHull pred:   ", pred_2[bad_ind])
        print("","%.1f%% difference between guesses."%(100.0*perc_diff))
        print()
        # Creating 3D visual of the surface over the data
        if test_points.shape[1] <= 2:
            # Produce a 3D/2D plot showing the difference in predictions
            p = Plot()
            p.add("Train Points",*(train_points.T),train_values,group=1)
            p.add("Test Point", [test_points[bad_ind][0]],
                  [test_points[bad_ind][1]], [pred_1[bad_ind]],group=2)
            p.add("Guess simplex", *(train_points[simp1].T),
                  train_values[simp1], shade=False,group=4)
            p.add("True simplex", *(train_points[simp2].T),
                  train_values[simp2], shade=False,group=5)
            p.add("Graphics Delaunay", *(train_points.T),
                  train_values, plot_type='surface', use_gradient=True,group=3)

            # Produce a 2D plot of the two simplices
            if test_points.shape[1] > 1:
                p2 = Plot()
                p2.add("Train Points",*(train_points.T),
                       group=1,show_in_legend=False)
                p2.add("Test Point", [test_points[bad_ind][0]],
                       [test_points[bad_ind][1]], group=2,
                       show_in_legend=False)
                p2.add("Guess simplex", *(train_points[simp1].T),
                       shade=False, group=4, show_in_legend=False,
                       mode='lines', fill='toself', opacity=0.2)
                p2.add("True simplex", *(train_points[simp2].T),
                       shade=False, group=5, show_in_legend=False,
                       mode='lines', fill='toself', opacity=0.2)
                # plot both together
                multiplot([[p2,p]])
            else:
                # Just plot the 3D
                p.plot()

    return perc_diff

if __name__ == "__main__":

    #      Read in the performance data files     
    # ============================================
    qhull_file = "Output-qHullDelaunay-MDA_results.csv"
    # qhull_file = "Output-NewDelaunay-MDA_results.csv"
    qhull_data = read_struct(qhull_file, verbose=True)
    print(qhull_data.shape)

    delaunay_file = "Output-FinalDelaunay-MDA_results.csv"
    delaunay_data = read_struct(delaunay_file, verbose=True)
    print(delaunay_data.shape)

    header = list(delaunay_data.dtype.names)

    tolerance = 0.01 # minimum difference allowed between qhull and delaunay
    diffs = [] # List of (% diff, q_row, d_row)
    target = "Mean_Guess"
    to_show = 3 # Only show this many 

    # Only collect comparison in 1 and 2 dimensions
    only_plottable = False
    to_test = -1
    to_skip = 0 # 20

    # Identify the differences between the algorithms' performance
    for i in range(len(qhull_data)):
        q_row = qhull_data[i]
        d_row = delaunay_data[i]
        q_pt = np.array([q_row[h] for h in header[:4]])
        d_pt = np.array([d_row[h] for h in header[:4]])
        if not np.all(q_pt == d_pt):
            print(i)
            print(q_row)
            print(d_row)
            exit()
        difference = abs(q_row[target] - d_row[target])
        avg = (q_row[target] + d_row[target]) / 2
        if (difference / avg) > tolerance:
            diffs.append( (difference / avg, q_row, d_row) )

    # Sort by the magnitude of the percent difference in guesses
    diffs.sort(key=lambda i: -i[0])

    # Show off the most different cases
    print()
    for (perc_diff,q_row,d_row) in diffs:
        if to_show <= 0: break
        difference = abs(q_row[target] - d_row[target])
        avg = (q_row[target] + d_row[target]) / 2
        print()
        print("QHull:     ",q_row[target], q_row)
        print("Delaunay:  ",d_row[target], d_row)
        print("Difference:",difference)
        print("Average:   ",avg)
        print("Perc_Diff: ",perc_diff)
        to_show -= 1
    print()

    # Reduce the set of differences
    if only_plottable:
        to_remove = []
        # Get the worst difference that exists under desired conditions
        for i,(perc_diff, q_row, d_row) in enumerate(diffs):
            defined = [q_row[h] for h in header[:4]].count(-1)
            if defined < 2: to_remove.append( i )
        for i in to_remove[::-1]: diffs.pop(i)

    # Cycle through the differences (starting with the worst)
    for (perc_diff, q_row, d_row) in diffs[:to_test]:
        print("="*50)
        # Get the dimension of the test data
        col_names = [h for h in header[:4] if q_row[h] != -1]
        num_dims = len(col_names)
        train_perc = str(q_row["Train_Percentage"])
        dir_name = "%i-"%num_dims + "--".join(col_names)
        dir_name = "/" + os.path.join("Users", "thomaslux", "Desktop",
                                      "VarSys-2016-Mean_MDA_Data_0", dir_name)
        # Search the folder (with different training %s) for the desired value
        for folder_name in os.listdir(dir_name):
            training_percentage = os.path.basename(folder_name).split("--")[0]
            if train_perc in training_percentage:
                dir_name = os.path.join(dir_name, folder_name)
                break

        # Now we have the directory that averaged to produce the difference point
        print()
        print("","Directory with difference %.1f%%:"%(100.0*perc_diff))
        print("  ",dir_name)

        # Get the array version of the point that was worst (un-normalized)
        max_diff_point = np.array([q_row[h] for h in col_names])
        print("","Point with the largest difference:")
        print("  ",max_diff_point, perc_diff)

        #      Find out which file actually has the bad spot     
        # =======================================================
        worst = -float('inf')
        worst_pair = (None,None)

        # Cycle through known directory, only look for training files
        all_files = [f for f in os.listdir(dir_name) if "Train" in f]
        print()
        for f_name in all_files:
            print("\r  %i of %i"%(all_files.index(f_name)+1, len(all_files)),end="")
            train = os.path.join(dir_name,f_name)
            test = train.replace("Train", "Test")
            difference = test_performance(train, test)
            if difference > worst:
                worst = difference
                worst_pair = (train, test)
        print()

        # Only both plotting if a point was found
        if worst_pair[0] != None:
            print("",worst)
            print("  ",worst_pair[0])
            print("  ",worst_pair[1])
            # Now plot the worst of the differences
            test_performance(*worst_pair, plot=True)

        print("="*50)
