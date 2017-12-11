# See why and where the prediction error is so large
#  - check 4D data first, it's bad in some spots
#  - verify the convex hull?
#  - look at plots of th 1D predictions, sanity check

import numpy as np
from util.plotly import *
from py_settings import *

# Preparing the data
normalize_points = True
random_points = True
random_perturbation = False
tight_plot = False
plot_results = True
autoshow_plot = True
perform_test = False

perturbation_multiplier = 1.0
plot_points = 4000
num_points = 10
dim = 2
multiplier = 10
rand_seed = 6 # 2

# from sklearn.linear_model import RANSACRegressor as sklearn_alg
# class Alg(sklearn_alg):
#     def __call__(self,x):
#         return self.predict(x)

algorithms = []
# algorithms += [ Alg() ]
algorithms += [ SVR() ] 
algorithms += [ MARS() ]
algorithms += [ LSHEP() ]
algorithms += [ Delaunay() ]
# algorithms += [ FitBoxMesh() ]
algorithms += [ MaxBoxMesh() ]
# algorithms += [ LinearModel() ]
algorithms += [ VoronoiMesh() ]
algorithms += [ MLPRegressor() ] 
# algorithms += [ NearestPlane() ]
algorithms += [ NearestNeighbor() ]

np.random.seed(rand_seed)

# Generate testing points in a grid pattern
plot_range = [[-0.3,1.3]]*dim # [[0,1]]*dim # 
if random_points:
    # Generate testing points randomly spaced
    points = np.random.random(size=(num_points**dim,dim)) * num_points
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

if (normalize_points and (len(points) > 1)):
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
if plot_results:
    p = Plot()
    # Plotting the results
    print()
    print("Plotting..")
    print()
    p.add("Raw data",*(points.T))

for s in algorithms:
    name = (str(s.__class__).split("'")[1]).split(".")[-1]
    start = time.time()
    s.fit(points[:,:-1].copy(), points[:,-1].copy(), *EXTRA_ARGS.get(name,[]))
    total = time.time()-start
    print("%s:"%name)
    print(" ","%.2e second fit"%total)

    if plot_results:
        use_gradient = len(algorithms) <= 1
        surf_x = np.meshgrid(*[
            np.linspace(lower,upper,plot_points**(1/dim))
            for (lower,upper) in plot_range])
        surf_x = np.array([x.flatten() for x in surf_x]).T
        plot_kwargs = {"marker_size":3, "marker_line_width":1,
                       "opacity":0.8, "plot_type":"surface",
                       "use_gradient":use_gradient, }
        if len(algorithms) <= 1: plot_kwargs.update({"opacity":1.0})
        if len(algorithms) >= 3: plot_kwargs.update({"plot_type":"scatter3d",
                                                     "marker_size":2})
        start = time.time()
        surf_y = s(surf_x.copy())
        total = time.time() - start
        print(" ","%.2e second prediction at %i points"%(total,len(surf_y)))
        print(" "," ","%.2e second per point"%(total/len(surf_y)))
        p.add(name, *(surf_x.T), surf_y,  **plot_kwargs, )
        print()

if plot_results: p.plot(show=autoshow_plot)
