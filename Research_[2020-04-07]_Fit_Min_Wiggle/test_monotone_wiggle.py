import numpy as np
from fraction import Fraction
from monotone import monotone_quintic_spline
from util.plot import Plot
from util.data import Data

# Compute the L2 of a function over a range.
def L2(f, low, upp):
    f2_int = (f.derivative(2)**2).integral()
    return f2_int(upp) - f2_int(low)

# Construct multiple different fits over x and y values
def compare_fits(x_vals, y_vals, plot=False):
    x_vals = list(map(Fraction, x_vals))
    y_vals = list(map(Fraction, y_vals))

    # Execute with different operational modes, (tests happen internally).
    functions = dict(
        quintic_interp = monotone_quintic_spline(x_vals, y_vals, method="poly"),
        quadratic_regres = monotone_quintic_spline(x_vals, y_vals, method="quad"),
        quadratic_fill = monotone_quintic_spline(x_vals, y_vals, method="fill", mids=2, ends=2),
    )

    # Compute the L2 of the second derivatives.
    wiggles = {name+"_L2":L2(functions[name],min(x_vals),max(x_vals))
               for name in functions}
    # Compute the second derivative values. and percentiles of them.
    ddf_pos = np.linspace(float(x_vals[0]), float(x_vals[-1]), num=100)
    for name in functions:
        ddf_vals = abs(np.array(functions[name].derivative(2)(ddf_pos), dtype=float))
        for p in (0, 25, 50, 75, 100):
            wiggles[name+f"_percentile_{p:03d}"] = np.percentile(ddf_vals,p)

    # Plot the data if the user wanted that.
    if plot:
        p = Plot()
        p.add("Points", x_vals, y_vals)
        # Add the functions
        for name in sorted(functions):
            p.add_func(name, functions[name], [min(x_vals), max(x_vals)])
        # Add the L2 functions.
        for i,name in enumerate(sorted(functions)):
            f = functions[name]
            f2_int = (f.derivative(2)**2).integral()
            p.add_func(name, f2_int, [min(x_vals), max(x_vals)],
                       color=i+1, dash="dot", group=0)
        p.show()

    # Return the results.
    return wiggles


# Run the tests.
def run_tests(seed=0, trials=10, Ns=(4,8,16,32,64)):
    np.random.seed(seed)
    all_data = Data()
    for N in Ns:
        for t in range(trials):
            x = np.random.random(size=(N,))
            y = np.random.random(size=(N,))
            x.sort(); y.sort()
            out = compare_fits(x,y)
            names = sorted(out)
            results = [N, t] + list(map(float,(out[n] for n in names)))
            all_data.append(results)
            if t == 0:
                all_data.names = ["Num nodes", "Trial number"] + \
                                 [n.title().replace("_"," ") for n in names]
                print(all_data)
                print(names)
            print(results,",")
        # Save a place-holder file for the Data object.
        all_data.save(f"seed-{seed}.csv")
        # Save some data for just this "N" in a CSV.
        np_data = np.array(all_data)
        np.savetxt(f"N{N}_trials-{trials}_seed-{seed}.csv",
                   np_data, delimiter=",", header=",".join(names), comments='')



run_tests(seed=0)



# N = 2**4
# # Show some examples.
# x = np.random.random(size=(N,))
# y = np.random.random(size=(N,))
# x.sort(); y.sort()
# out = compare_fits(x, y, plot=True)
# 
    # from util.plot import Plot
    # plot_range = [min(x_vals)-.1, max(x_vals)+.1]
    # p = Plot()
    # p.add("Points", list(map(float,x_vals)), list(map(float,y_vals)))
    # p.add_func("f (min_d1=0)", f, plot_range)
    # p.add_func("f' (min_d1=0)", f.derivative(1), plot_range, dash="dash")
    # p.add_func("f'' (min_d1=0)", f.derivative(2), plot_range, dash="dot")
    #     print("L2(f''):", float(L2(f.derivative(2), x_vals[0], x_vals[-1])))
# 
# p = Plot()


# from util.stats import rank_probability
# 
# for i in range(len(names)):
#     print()
#     print(f"{names[i]}: ")
#     vals = [
#         np.percentile(all_data[:,i],v) for v in (0,25,50,75,100)
#     ]
#     print(f"  "+("%.3f  "*len(vals))%tuple(vals))
#     this_column = all_data[:,i]
#     other_columns = [all_data[:,j] for j in range(all_data.shape[1]) if j != i]
#     print(" is best %.0f%% of the time"%(100*rank_probability(
#         this_column, other_columns, order=True)))


# N = 6
# trials = 100
# 
#   min    25%   median   75%    max
# 
# f_(quadratic interpolant): 
#   0.028  0.133  0.205  0.296  0.588  
#  is best 22.0% of the time
# 
# f_(use a quintic): 
#   0.028  0.136  0.199  0.297  0.588  
#  is best 42.0% of the time
# 
# f_(weighted quadratic regression): 
#   0.026  0.132  0.201  0.291  0.587  
#  is best 36.0% of the time
