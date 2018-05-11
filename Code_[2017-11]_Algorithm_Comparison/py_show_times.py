from util.plotly import *
import pickle
from py_settings import *

# Get the names of the series, order them for aesthetics
with open(FILE_NAME,"rb") as f:
    fit_times, eval_times = pickle.load(f)
names = list(fit_times.keys())

# # Fixing bad data
# index = 0
# while (fit_times[names[0]][index][0] <= 9310): index += 1
# for n in names:
#     fit_times[n] = fit_times[n][:index]
#     eval_times[n] = eval_times[n][:index]
# with open(FILE_NAME,"wb") as f:
#     pickle.dump((fit_times,eval_times),f)
# exit()

# Sort by (longest length first, longest execution time first)
names.sort( key=lambda n: (-len(fit_times[n]),-(eval_times[n][-1][1]+fit_times[n][-1][1])) )

# Setup title based on measurements
if MEASURE_PTS:
    first = "increasing numbers of points"
    second = "%i dimensions"%DIM
    x_title = "Number of points"
elif MEASURE_DIM:
    first = "%i points"%NUM_POINTS
    second = "increasing dimension"
    x_title = "Number of dimensions"
title = "Fitting %s in %s, evaluating at %i points"%(
    first,second,NUM_EVALS)



p = Plot(title, x_title, "Execution Time (seconds)")

for n in names:
    p.color_num += 1
    alpha = 0.5
    eval_brightness = 1.0
    fit_brightness = 0.7

    color = p.color(p.color_num, alpha=alpha, brightness=eval_brightness)
    p.add(n+" Eval", *list(zip(*eval_times[n])), mode='lines',
          fill='tonexty', color=color, group=n)
    color = p.color(p.color_num, alpha=alpha, brightness=fit_brightness)
    p.add(n+" Fit", *list(zip(*fit_times[n])), mode='lines',
          fill='tozeroy', color=color, group=n)

    if n == "Delaunay":
        from scipy.optimize import minimize
        print("FITTING")
        # Collect the dimension (x) and total time (y) into arrays
        dims = np.array(fit_times[n])
        dims, times = dims[:,0], dims[:,1]
        times += np.array(eval_times[n])[:,1]
        print(dims)
        print(times)
        coefs = np.zeros(4)
        objective = lambda coef: np.sum((np.polyval(coef, dims) - times)**2)
        solution = minimize(objective, coefs)
        print(solution)        
        solution = solution.x
        p.add_func(n+" fit", lambda x: np.polyval(solution, x), [min(dims),max(dims)*3])



p.plot()
