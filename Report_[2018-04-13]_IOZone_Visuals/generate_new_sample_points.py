from util.data import Struct

data = Struct().load("prediction_results.csv")
data.summarize(15)

config_names = data.names[:-1]
uniques = {n:sorted(set(v for v in data[n])) for n in config_names}
for n in uniques:
    print(n, sorted(uniques[n]))

# Initialize a structure for holding all of the data
box_errors = Struct(names=["Expected Midpoint", "Actual Midpoint", "Num Sources", "Error"])

import os
if (not os.path.exists("grid_box_errors.pkl.gz")):
    for freq_i in range(len(uniques["Frequency"])-1):
        print(freq_i+1,":",len(uniques["Frequency"])-1)
        for file_i in range(len(uniques["File Size"])-1):
            print("  ", file_i+1,":", len(uniques["File Size"])-1)
            for rec_i in range(len(uniques["Record Size"])-1):
                for thread_i in range(len(uniques["Num Threads"])-1):
                    total_error = 0
                    freqs = {uniques["Frequency"][freq_i], uniques["Frequency"][freq_i+1]}
                    files = {uniques["File Size"][file_i], uniques["File Size"][file_i+1]}
                    recs  = {uniques["Record Size"][rec_i], uniques["Record Size"][rec_i+1]}
                    thrds = {uniques["Num Threads"][thread_i], uniques["Num Threads"][thread_i+1]}
                    square_data = data[data["Frequency"] == freqs]
                    square_data = square_data[square_data["File Size"] == files]
                    square_data = square_data[square_data["Record Size"] == recs]
                    square_data = square_data[square_data["Num Threads"] == thrds]
                    total_error = sum(square_data["KS Statistic for Prediction"])
                    correct_midpoint = (
                        sum(freqs) / 2,
                        sum(files) / 2,
                        sum(recs)  / 2,
                        sum(thrds) / 2
                    )
                    actual_midpoint = tuple()
                    if len(square_data) > 0:
                        actual_midpoint = (sum(square_data["Frequency"]) / len(square_data),
                                           sum(square_data["File Size"]) / len(square_data),
                                           sum(square_data["Record Size"]) / len(square_data),
                                           sum(square_data["Num Threads"]) / len(square_data),
                        )
                    box_errors.append([correct_midpoint, actual_midpoint, 
                                       len(square_data), total_error])

    # Save the data that was generated
    box_errors.save("grid_box_errors.pkl.gz")
else:
    box_errors.load("grid_box_errors.pkl.gz")

print()

# Summarize the data (excluding the actual configurations)
box_errors.summarize(15)
                
# print(box_errors[ box_errors[box_errors.names[-1]] > .4 ])
# print(len(set(box_errors["Config"])))


# Compilation started at Tue Apr 17 11:22:39
# 
# python3 generate_new_sample_points.py 
# SUMMARY:
# 
#   This data has 4052 rows, 5 columns.
#     5 columns are recognized as ordinal, 0 are categorical.
#     0 rows have missing values.
# 
# COLUMNS:
# 
#   0 -- "Frequency"                   <class 'int'>   (7 unique values)
#     1200000  531 ( 13.1%) #######
#     1600000  531 ( 13.1%) #######
#     2000000  531 ( 13.1%) #######
#     2300000  605 ( 14.9%) #######
#     2800000  756 ( 18.7%) #########
#     3200000  549 ( 13.5%) #######
#     3500000  549 ( 13.5%) #######
# 
#   1 -- "File Size"                   <class 'int'>   (11 unique values)
#     4         63 (  1.6%) #
#     16       162 (  4.0%) ##
#     64       261 (  6.4%) ###
#     256      360 (  8.9%) ####
#     1024     443 ( 10.9%) #####
#     4096     513 ( 12.7%) ######
#     8192     576 ( 14.2%) #######
#     16384    639 ( 15.8%) ########
#     32768    360 (  8.9%) ####
#     65536    360 (  8.9%) ####
#     1048576  315 (  7.8%) ####
# 
#   2 -- "Record Size"                 <class 'int'>   (13 unique values)
#     4      639 ( 15.8%) ########
#     8      225 (  5.6%) ###
#     16     576 ( 14.2%) #######
#     32     189 (  4.7%) ##
#     64     513 ( 12.7%) ######
#     128    146 (  3.6%) ##
#     256    450 ( 11.1%) ######
#     512    108 (  2.7%) #
#     1024   387 (  9.6%) #####
#     2048    81 (  2.0%) #
#     4096   324 (  8.0%) ####
#     8192   216 (  5.3%) ###
#     16384  198 (  4.9%) ##
# 
#   3 -- "Num Threads"                 <class 'int'>   (9 unique values)
#     1     451 ( 11.1%) ######
#     8     451 ( 11.1%) ######
#     16    450 ( 11.1%) ######
#     24    450 ( 11.1%) ######
#     32    450 ( 11.1%) ######
#     40    450 ( 11.1%) ######
#     48    450 ( 11.1%) ######
#     56    450 ( 11.1%) ######
#     64    450 ( 11.1%) ######
# 
#   4 -- "KS Statistic for Prediction" <class 'float'> (3861 unique values)
#     [0.00e+00, 7.14e-02)  448 ( 11.1%) ######
#     [7.14e-02, 1.43e-01)  713 ( 17.6%) #########
#     [1.43e-01, 2.14e-01)  426 ( 10.5%) #####
#     [2.14e-01, 2.86e-01)  484 ( 11.9%) ######
#     [2.86e-01, 3.57e-01)  553 ( 13.6%) #######
#     [3.57e-01, 4.29e-01)  638 ( 15.7%) ########
#     [4.29e-01, 5.00e-01)  228 (  5.6%) ###
#     [5.00e-01, 5.71e-01)  184 (  4.5%) ##
#     [5.71e-01, 6.43e-01)  184 (  4.5%) ##
#     [6.43e-01, 7.14e-01)  133 (  3.3%) ##
#     [7.14e-01, 7.86e-01)   38 (  0.9%) 
#     [7.86e-01, 8.57e-01)   24 (  0.6%) 
#     [8.57e-01, 9.29e-01)    8 (  0.2%) 
#     [9.29e-01, 1.00e+00]   15 (  0.4%) 
# 
# Frequency [1200000, 1600000, 2000000, 2300000, 2800000, 3200000, 3500000]
# File Size [4, 16, 64, 256, 1024, 4096, 8192, 16384, 32768, 65536, 1048576]
# Record Size [4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384]
# Num Threads [1, 8, 16, 24, 32, 40, 48, 56, 64]
# 
# SUMMARY:
# 
#   This data has 5760 rows, 4 columns.
#     2 columns are recognized as ordinal, 2 are categorical.
#     0 rows have missing values.
# 
# COLUMNS:
# 
#   0 -- "Expected Midpoint" <class 'tuple'> (5760 unique values)
#     (1400000.0, 10.0, 6.0, 4.5)             1 (  0.0%) 
#     (1400000.0, 10.0, 6.0, 12.0)            1 (  0.0%) 
#     (1400000.0, 10.0, 6.0, 20.0)            1 (  0.0%) 
#     (1400000.0, 10.0, 6.0, 28.0)            1 (  0.0%) 
#     (1400000.0, 10.0, 6.0, 36.0)            1 (  0.0%) 
#     (1400000.0, 10.0, 6.0, 44.0)            1 (  0.0%) 
#     (1400000.0, 10.0, 6.0, 52.0)            1 (  0.0%) 
#     (1400000.0, 10.0, 6.0, 60.0)            1 (  0.0%) 
#     (1400000.0, 10.0, 12.0, 4.5)            1 (  0.0%) 
#     (1400000.0, 10.0, 12.0, 12.0)           1 (  0.0%) 
#     (1400000.0, 10.0, 12.0, 20.0)           1 (  0.0%) 
#     (1400000.0, 10.0, 12.0, 28.0)           1 (  0.0%) 
#     (1400000.0, 10.0, 12.0, 36.0)           1 (  0.0%) 
#     (1400000.0, 10.0, 12.0, 44.0)           1 (  0.0%) 
#     (1400000.0, 10.0, 12.0, 52.0)           1 (  0.0%) 
#     ... (increase 'max_display' to see more summary statistics).
# 
#   1 -- "Actual Midpoint"   <class 'tuple'> (3491 unique values)
#     ()                                                                               1392 ( 24.2%) ############
#     (1400000.0, 16.0, 16.0, 4.5)                                                        2 (  0.0%) 
#     (1400000.0, 16.0, 16.0, 12.0)                                                       2 (  0.0%) 
#     (1400000.0, 16.0, 16.0, 20.0)                                                       2 (  0.0%) 
#     (1400000.0, 16.0, 16.0, 28.0)                                                       2 (  0.0%) 
#     (1400000.0, 16.0, 16.0, 36.0)                                                       2 (  0.0%) 
#     (1400000.0, 16.0, 16.0, 44.0)                                                       2 (  0.0%) 
#     (1400000.0, 16.0, 16.0, 52.0)                                                       2 (  0.0%) 
#     (1400000.0, 16.0, 16.0, 60.0)                                                       2 (  0.0%) 
#     (1400000.0, 40.0, 16.0, 4.5)                                                        2 (  0.0%) 
#     (1400000.0, 40.0, 16.0, 12.0)                                                       2 (  0.0%) 
#     (1400000.0, 40.0, 16.0, 20.0)                                                       2 (  0.0%) 
#     (1400000.0, 40.0, 16.0, 28.0)                                                       2 (  0.0%) 
#     (1400000.0, 40.0, 16.0, 36.0)                                                       2 (  0.0%) 
#     (1400000.0, 40.0, 16.0, 44.0)                                                       2 (  0.0%) 
#     ... (increase 'max_display' to see more summary statistics).
# 
#   2 -- "Num Sources"       <class 'int'>   (12 unique values)
#     0    1392 ( 24.2%) ############
#     4     544 (  9.4%) #####
#     6      56 (  1.0%) 
#     8    2020 ( 35.1%) ##################
#     9       2 (  0.0%) 
#     10    246 (  4.3%) ##
#     11      2 (  0.0%) 
#     12    550 (  9.5%) #####
#     13      2 (  0.0%) 
#     14     46 (  0.8%) 
#     15      2 (  0.0%) 
#     16    898 ( 15.6%) ########
# 
#   3 -- "Error"             <class 'float'> (3491 unique values)
#     [0.00e+00, 7.32e-01) 1910 ( 33.2%) #################
#     [7.32e-01, 1.46e+00)  709 ( 12.3%) ######
#     [1.46e+00, 2.20e+00)  694 ( 12.0%) ######
#     [2.20e+00, 2.93e+00)  665 ( 11.5%) ######
#     [2.93e+00, 3.66e+00)  527 (  9.1%) #####
#     [3.66e+00, 4.39e+00)  399 (  6.9%) ###
#     [4.39e+00, 5.13e+00)  262 (  4.5%) ##
#     [5.13e+00, 5.86e+00)  173 (  3.0%) ##
#     [5.86e+00, 6.59e+00)  161 (  2.8%) #
#     [6.59e+00, 7.32e+00)   88 (  1.5%) #
#     [7.32e+00, 8.05e+00)   72 (  1.2%) #
#     [8.05e+00, 8.79e+00)   55 (  1.0%) 
#     [8.79e+00, 9.52e+00)   31 (  0.5%) 
#     [9.52e+00, 1.03e+01]   13 (  0.2%) 
# 
# 
# Compilation finished at Tue Apr 17 11:22:40
