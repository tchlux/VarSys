import numpy as np
from fits import cdf_points, flat_fit, linear_fit, cubic_fit, quintic_fit
from util.plot import Plot
from random import seed, sample

values = list(np.percentile(np.random.normal(size=(100000,)),
                            np.linspace(0,100,1000)))
truth = linear_fit(values)
true_min_max = (min(values), max(values))

# Create a visual of some sample distribution approximations.
def make_plot(functions, prename):
    # Initialize some settings.
    seed(0); k = 10
    pop = sample(values, k)
    x, y = cdf_points(pop)
    styles = [None, "dashdot", "dot", "dash"]
    styles = [None] * 4
    # Create the plot.
    p = Plot("", "x","CDF", font_family="times", font_size=18)
    p.add("Sample", x, y)
    p.add_func("Truth", truth, true_min_max, color=p.color((0,0,0,.3)))
    for f,s in zip(functions, styles):
        name = f.__name__.replace("_"," ").title().split()[0].replace("Flat","EDF")
        # Set the legend properties.
        if "quintic" in name.lower():
            p.add_func(name, f(pop), true_min_max, dash=s, opacity=.8, fill='toprevy')
        else:
            p.add_func(name, f(pop), true_min_max, dash=s, opacity=.8)
    legend = dict(
        xanchor = "center",
        yanchor = "top",
        x = .25,
        y = .8,
        orientation = "v",
        bgcolor="white",
        bordercolor="grey",
        borderwidth=.5
    )
    # Create the plot.
    # p.show(y_range=[-.1, 1.1], x_range=true_min_max, width=400*1.4,
    #        height=300*1.4) #, file_name=prename+"-sample-prediction.html")

    # - remove the stuff from quintic fit
    # - remove error bar plots


    # Fit the errors
    p = Plot("", "Absolute Error","CDF", font_family="times", font_size=18)
    print("Computing errors..")
    fit = f(pop)
    errors = abs(np.linspace(0,1,len(values)) - np.array([fit(v) for v in values]))
    print("Fitting error distribution..")
    fit = linear_fit(errors)
    print("Making plot..")
    p.add_func("Error", fit, [min(errors),max(errors)])
    p.show(width=400*1.4, height=300*1.4)

    

# Make the two different plots.
functions = [flat_fit, linear_fit, cubic_fit, quintic_fit]
make_plot(functions, "fl")

# functions = [cubic_fit, quintic_fit]
# make_plot(functions, "cq")
