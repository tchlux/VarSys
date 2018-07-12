import autograd.numpy as np
from util.plot import Plot

legend = dict(
    xanchor = "center",
    yanchor = "top",
    x = .5,
    y = 1.1,
    orientation = "h",
)

f = lambda x: np.sin(x**(1/2))

# ====================================================================
# Neural network training demonstration
from util.algorithms.neural_network import MLPRegressor as MLP

# Calculate training points
num_pts = 7
x = np.linspace(1,(4*np.pi)**2,num_pts)
y = f(x).flatten()

# Construct model
hidden_nodes = 100
model = MLP(hidden_layer_sizes=(hidden_nodes,), solver='sgd',
            max_iter=3000, random_state=10, tol=-float('inf'),
            verbose=True, shuffle=False)
model.fit(x[:,None] , y)
predict = lambda x: model.predict(x)

# Make plot of Neural Network fit
p = Plot("",font_size=20)
p.add_func("Truth Function", f, [min(x), max(x)], dash="dash",
           color=p.color(1, alpha=.25))
p.add_func("Neural Network", predict, [min(x), max(x)],
           color=p.color(0), vectorized=True)
p.add("Training Points", x, y, color=p.color(7))
p.show(file_name="MLP_fit.html", legend=legend)


# --------------------------------------------------------------------
# Delaunay interpolation demonstration
from util.algorithms import qHullDelaunay as Delaunay

num_pts = 7
x = np.linspace(1,(4*np.pi)**2,num_pts).reshape((num_pts,1))
y = f(x).flatten()

model = Delaunay()
model.fit(x, y)

# Make plot of Delaunay fit  f"Delaunay Fit of {num_pts} Points"
p = Plot("",font_size=20)
p.add_func("Truth Function", f, [min(x), max(x)], dash="dash",
           color=p.color(1, alpha=.25))
p.add_func("Delaunay", model, [min(x)-10, max(x)+10], color=p.color(0))
p.add("Training Points", x[:,0], y, color=p.color(7))
p.show(file_name="Delaunay_fit.html", legend=legend)

# ====================================================================
