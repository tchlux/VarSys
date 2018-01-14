from util.plotly import Plot, multiplot
import pickle, os
from util.stats import normal_pdf, normal_cdf, cdf_fit_func, plot_percentiles

# Plot a normal PDF and CDF
p1 = Plot("", "X value", "P[x = X]")
p1.add_func("PDF of N(0,1)", normal_pdf, [-4,4])
p2 = Plot("Probability Density Function (PDF) and Cumulative Distribution Function (CDF) <br> of Normal Distribution with zero mean and unit variance",
          "X value", "P[x < X]")
p2.add_func("CDF of N(0,1)", normal_cdf, [-4,4], color=p2.color(1))
multiplot([[p1],[p2]], shared_x=True, file_name="Normal_Distribution_PDF_CDF.html")


# Generate data (approximation of normal distribuitons with increasing samples)
if not os.path.exists("plot_info.pkl"):
    p = Plot()
    num_normals = 100
    percentiles = np.linspace(0,100,11)
    x = np.random.normal(size=(num_normals,10))
    for size in np.linspace((10)**(1/2),(1000)**(1/2),100):
        size = int(round(size**2))
        if (size > len(x[0])):
            needed = size - len(x[0])
            x = np.concatenate((x, np.random.normal(size=(num_normals,needed))), axis=1)
        print("Starting",size)
        funcs = [cdf_fit_func(row) for row in x]
        test_x = np.linspace(-4,4,100)
        test_y = np.array([f(test_x) for f in funcs]).T
        plot_percentiles(p, "Normal CDF", test_x, test_y, percentiles, frame=size)
    with open("plot_info.pkl", "wb") as f:
        pickle.dump(p, f)
else:
    with open("plot_info.pkl", "rb") as f:
        p = pickle.load(f)

p.title = "Percentile clouds of possible CDFs from 100 normal distributions with N samples"
p.y_title = "Range of empirically estimated CDFs"
p.x_title = "X value (0 mean, 1 stdev)"
p.plot(frame_prefix="Number of Samples: ", show_legend=False, 
       file_name="Convergence_of_CDFs_Example.html")


