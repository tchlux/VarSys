
import numpy as np
from util.misc.grid import grid as grid_func
from util.plot import Plot
from util.math import transpose

np.random.seed(0)
# x = np.random.random((3,2))

names = ["equilateral", "right", "skinny"]
xs = [
    np.array([
        [.45,.7],
        [.6, .3],
        [.3,.3],
    ]),
    np.array([
        [.45, .7],
        [.45,.3],
        [.3,.3],
    ]),
    np.array([
        [.45, .7],
        [.45,.3],
        [.4,.5],
    ])
]

for (name, x) in zip(names, xs):

    space = .2
    grid = grid_func(np.min(x)-space, np.max(x)+space, dim=2, num=10000)
    # grid += np.random.random(grid.shape) * .05
    # grid = np.random.random((100,2))
    ones = np.ones((x.shape[0],1))
    labels = []

    # print("x.shape: ",x.shape)
    # print("ones:    \n",ones)

    for pt in grid:
        system = np.concatenate((x, ones), axis=1).T
        x_pt = np.concatenate((pt,[1]))
        # print("system:  ",system)
        # print("x_pt:    ",x_pt)
        weights = np.linalg.solve(system, x_pt)
        # print("weights: ",weights, sum(weights))
        labels.append( "  ".join([f"{w:.2f}" for w in weights]) )
        # print(labels[-1])
        # exit()

    p = Plot(f"Affine weights for {name} triangle")
    p.add("grid", *(grid.T), text=labels, color=p.color(color=(255,255,255)),
          opacity=0, hoverinfo="text")
    p.add("Points", *transpose(list(x) + [x[0]]), marker_size=20,
          mode="markers+lines+text", text=["1", "2", "3", ""],
          line_color="rgba(0,0,0,.5)", shade=False)
    p.graph(hovermode="closest", file_name="affine_weights.html",
            show=(name==names[-1]), append=(name!=names[0]))
