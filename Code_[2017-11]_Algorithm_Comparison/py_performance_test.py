import numpy as np
from util.plotly import *
import pickle

from py_settings import *

fit_times = {}
eval_times = {}
max_iters = MAX_ITERS
surfs = [NearestNeighbor(), LinearModel(), FitBoxMesh(),
         qHullDelaunay(), Delaunay(), LSHEP(), MARS(),
         BayesTree(), NearestPlane(), MaxBoxMesh(), tgp(),
         dynaTree(), BBS(), MLPRegressor()]
num_points = NUM_POINTS
dim = DIM
surfs = surfs[:4]
max_iters = 100

if __name__ == "__main__":
    while (max_iters > 0) and (len(surfs) > 1):
        max_iters -= 1
        if MEASURE_PTS:
            num_points += 100
        elif MEASURE_DIM:
            dim += 1

        # Generate testing points randomly spaced
        points = np.random.random(size=(num_points,dim)) * (num_points-1)
        # Calculate the associated response values
        values = np.array([FUN(pt) for pt in points])

        print()
        print("Points:   %s"%(str(points.shape)))
        print("Response: %s"%(str(values.shape)))
        print("Beginning to fit the data...")

        surfs = [Delaunay(), MaxBoxMesh()]

        to_remove = []
        for s in surfs:
            name = GET_CLASS_NAME(s)
            # Test time-to-fit
            start = time.time()
            s.fit(points, values, *EXTRA_ARGS.get(name,[]))
            fit_time = time.time()-start
            # Test time-to-evaluate
            start = time.time()
            s( points[:NUM_EVALS] )
            eval_time = time.time()-start
            # Depending on what we're measuring, record appropriately
            idx = 0 if MEASURE_PTS else 1
            fit_times[name] = fit_times.get(name,[]) + [[points.shape[idx],fit_time]]
            eval_times[name] = eval_times.get(name,[]) + [[points.shape[idx],eval_time]]
            # Calculate the total time to give to the user
            total = fit_time + eval_time
            print("%s: %s seconds"%(name,total))
            # Mark algorithms that take more than "max_compute_time_sec"
            if total > MAX_COMPUTE_TIME_SEC:
                to_remove.append(s)

        # Remove algorithms that are taking too long to compute
        for s in to_remove:
            surfs.remove(s)

    with open(FILE_NAME,"wb") as f:
        pickle.dump((fit_times,eval_times),f)
