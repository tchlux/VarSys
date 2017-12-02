from pyasis.algorithms import Delaunay
from pyasis.plotly import Plot
import numpy as np

x = np.array([1,3,3.5,5,6,7])
y = np.array([1,0,1.8,2.5,5.5,6])
min_max = [min(x), max(x)]
shift = min_max[0]
scale = min_max[1] - min_max[0]

# Define a truth function
fun = lambda new_x: np.interp(new_x, x, y)

# Pick a number of samples
num_samples = 8

# Generate the evenly and logarithmically spaced sample points
even_x = np.linspace(min(x),max(x),num_samples)
even_y = fun(even_x)
even_x = even_x.reshape((num_samples, 1))
even_d = Delaunay()
even_d.fit(even_x, even_y)
even_d_out = lambda x: even_d(x)

# Fit the exponentially spaced points
log_x = np.linspace(np.log(min(x)), np.log(max(x)), num_samples)
log_y = np.log( fun(np.exp(log_x)) )
log_x = log_x.reshape((num_samples,1))
# Generate a non-trasnformed fit
log_d = Delaunay()
log_d.fit(np.exp(log_x), np.exp(log_y))
log_d_out = lambda x: log_d(x)
# Generate a transformed fit
log_d_trans = Delaunay()
log_d_trans.fit(log_x, log_y)
log_d_trans_out = lambda x: np.exp( log_d_trans( np.log(x) ) )

even_color = "rgb(30,180,80)"
log_color = "rgb(100,50,200)"

p = Plot("The Dangers Associated with Log-Scaling")
p.add_func("Truth", fun, min_max, group="Truth")
p.add("Evenly Spaced Points", even_x[:,0], even_y, color=even_color)
p.add_func("Delaunay Even", even_d_out, min_max, color=even_color)
p.add("Exponentially Spaced Points", np.exp(log_x[:,0]), np.exp(log_y), color=log_color)
p.add_func("Delaunay Log Normal", log_d_out, min_max, color=log_color)
p.add_func("Delaunay Log Transformed", log_d_trans_out, min_max, color=log_color)

p.plot()
