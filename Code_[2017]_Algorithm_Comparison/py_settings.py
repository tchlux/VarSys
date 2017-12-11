from util.decorators import timeout_after
from util.classes import AtomicOpen

import numpy as np
# For retrieving the comparing algorithms
import os, sys, time, fcntl, signal, inspect
from util.algorithms import *

# Get the name of a class as listed in python source file.
GET_CLASS_NAME = lambda obj: (repr(obj)[1:].split(" ")[0]).split(".")[-1]

# Reduce a number to its minimal display form (no unnecessary 0's on right)
CLEAN_NUMBER_STRING = lambda number: str(number).rstrip("0").rstrip(".") \
                      if "." in str(number) else str(number)
# Number used for dividing error (to get relative error) if true value is 0
SMALL_NUMBER = 0.00001
# Arguments for algorithms that need to be controlled
FBM_NUM_BOXES = MARS_MAX_BASES = 200
NEAREST_NEIGHBOR_NUM_NEIGHBORS = 1
BAYES_NUM_TREES = 100
EXTRA_ARGS = {
    "FitBoxMesh": [FBM_NUM_BOXES],
    "MARS": [MARS_MAX_BASES],
    "BayesTree": [BAYES_NUM_TREES],
    "NearestNeighbor": [NEAREST_NEIGHBOR_NUM_NEIGHBORS]
}
OUTPUT_DATA_FILE = "MDA_results.csv"

FIT_TIMEOUT_SECONDS = 60
APPROX_TIMEOUT_SECONDS = 60
FAIL_VALUE = "FAILED"

# Preparing the data
MEASURE_PTS   = False
MEASURE_DIM   = not MEASURE_PTS

# FUN = lambda x: 100 * np.sum(x)**2
# FUN = lambda x: (x[0]-2)*(x[1]-1) if ((x[0] > 2) and (x[1] > 1)) else 0
# FUN = lambda x: 1*x[0]
# FUN = lambda x: (x[0]-num_points/3)*(x[1]-num_points/2)**2
FUN = lambda x: np.cos(x[0]) * (np.sin(x[1]) if len(x) > 1 else 1.0) + 1
# FUN = lambda x: np.product(np.sin(np.array(x))/10)

MAX_COMPUTE_TIME_SEC = 5.0
MAX_ITERS = 30
NUM_EVALS = 10
PLOT_POINTS = 1000

# Settings determined based on what we're measuring
if MEASURE_PTS:
    NUM_POINTS = 10
    DIM = 5
    FILE_NAME = "Run-time-pts.pkl"
elif MEASURE_DIM:
    NUM_POINTS = 500
    DIM = 1
    FILE_NAME = "Run-time-dim.pkl"

