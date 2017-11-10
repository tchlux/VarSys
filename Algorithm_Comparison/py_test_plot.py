# See why and where the prediction error is so large
#  - check 4D data first, it's bad in some spots
#  - verify the convex hull?
#  - look at plots of th 1D predictions, sanity check

import numpy as np
from util.plotly import *
from py_settings import *

# Preparing the data
normalize_points = True
random_points = False
random_perturbation = False
tight_plot = False
plot_results = True
perform_test = False

perturbation_multiplier = 1.0
plot_points = 4000
num_points = 10
dim = 2
multiplier = 10

# from sklearn.linear_model import RANSACRegressor as sklearn_alg
# class Alg(sklearn_alg):
#     def __call__(self,x):
#         return self.predict(x)

# algorithms = [Alg()]
# algorithms = [BSB()]
# algorithms = [FitBoxMesh()]
# algorithms = [MaxBoxMesh()]
# algorithms = [SplineMesh()]
# algorithms = [NearestNeighbor()]
# algorithms = [LSHEP()]
# algorithms = [LinearModel()]
# algorithms = [Delaunay()]
# algorithms = [ MLPRegressor() ] 
algorithms = [ SVR() ] 

# algorithms = [NearestPlane()]

np.random.seed(0)

# Generate testing points in a grid pattern
plot_range = [[-0.3,1.3]]*dim # [[0,1]]*dim # 
if random_points:
    # Generate testing points randomly spaced
    points = np.random.random(size=(num_points**dim,dim)) * (num_points-1)*num_points
else:
    points = np.meshgrid(*[np.linspace(0,num_points-1,num_points) for i in range(dim)])
    points = np.array([p.flatten() for p in points]).T * num_points

if random_perturbation:
    points += np.random.random(size=points.shape)*perturbation_multiplier

FUN = lambda x: 100 * np.sum(x)**2
# FUN = lambda x: (x[0]-2)*(x[1]-1) if ((x[0] > 2) and (x[1] > 1)) else 0
# FUN = lambda x: 1*x[0]
# FUN = lambda x: (x[0]-num_points/3)*(x[1]-num_points/2)**2
# FUN = lambda x: np.cos(x[0]) * (np.sin(x[1]) if len(x) > 1 else 1.0) + 1
# FUN = lambda x: np.product(np.sin(np.array(x))/10)

points *= 0.05
# Calculate the associated response values
values = np.array([[FUN(pt) for pt in points]]).T
points = np.concatenate((points, values), axis=1)
# Sort the points so they resemble grid order
points = points[points[:,1].argsort()]
points = points[points[:,0].argsort()]

print()
print("Points:   %s"%(str(points[:,:-1].shape)))#,points[:10,:-1])
print("Response: %s"%(str(points[:, -1].shape)))#,points[:10,-1])

if normalize_points:
    # Normalize the points themselves
    max_val = float(np.max(points[:,:-1]))
    min_val = float(np.min(points[:,:-1]))
    points[:,:-1] = (points[:,:-1] - min_val) / (max_val - min_val)
    # Normalize the response values
    max_val = float(np.max(points[:,-1]))
    min_val = float(np.min(points[:,-1]))
    points[:,-1] = (points[:,-1] - min_val) / (max_val - min_val) + 2.0
else:
    min_val = np.min(points, axis=0)
    max_val = np.max(points, axis=0)
    plot_range = [[ plot_range[i][0] * (max_val[i]-min_val[i]) + min_val[i],
                    plot_range[i][1] * (max_val[i]-min_val[i]) + min_val[i] ]
                  for i in range(dim)]

if tight_plot:
    plot_range = [[min(points[:,i]),max(points[:,i])]
                  for i in range(dim)]


print("Beginning to fit the data...")
for s in algorithms:
    name = (str(s.__class__).split("'")[1]).split(".")[1]
    start = time.time()
    s.fit(points[:,:-1].copy(), points[:,-1].copy(), *EXTRA_ARGS.get(name,[]))
    total = time.time()-start
    print("%s: %s seconds"%(name,total))


# surf = algorithms[0]
# print(surf.boxes.T)
# print(surf.box_widths.T)
# exit()

if plot_results:
    # Plotting the results
    print()
    print("Plotting..")
    print()
    p = Plot()
    p.add("Raw data",*(points.T))
    use_gradient = len(algorithms) <= 1
    surf_x = np.meshgrid(*[
        np.linspace(lower,upper,PLOT_POINTS**(1/dim))
        for (lower,upper) in plot_range])
    surf_x = np.array([x.flatten() for x in surf_x]).T
    marker_args = {"marker_size":3, "marker_line_width":1, "opacity":0.8}
    if len(algorithms) <= 1: marker_args = {}
    for s in algorithms:
        name = (str(s.__class__).split("'")[1]).split(".")[1]
        print("Adding '%s'..."%(name))
        start = time.time()
        surf_y = s(surf_x.copy())
        print(surf_x.shape)
        print(surf_y.shape)
        total = time.time() - start
        print("  %.2f seconds."%total)
        p.add(name, *(surf_x.T), surf_y, use_gradient=use_gradient,
              plot_type='surface',
              # mode="markers", marker_size=3, marker_line_width=1,
              **marker_args,
        )
    p.plot()

if perform_test:
    print()
    print("Input points:")
    print(points[:,:-1])
    print()
    tester = Delaunay()
    tester.fit(points[:,:-1], points[:,-1])
    points = np.array([
        [0.5,1.0],
        [0.6, 1.2],
        [.6455161290322581902, 1.354483870967742032],
    ])

    tester.find_simplex(points)
