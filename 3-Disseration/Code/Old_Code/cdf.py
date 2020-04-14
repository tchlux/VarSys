if __name__ == "__main__":
    import numpy as np
    from util.plot import Plot
    from util.stats import cdf_points
    from monotone import monotone_quintic_spline

    np.random.seed(0)
    N = 100
    x = np.random.normal(size=(N,))

    points = np.array(list(zip(*cdf_points(x))))
    step = N // 3

    p = Plot()
    funcs = []
    p.color_num += 1
    for i in range(step):
        indices = np.arange(i,len(points),step)
        sub_points = points[indices]
        print()
        print(f"Constructing function {i+1}..", end=" ", flush=True)
        funcs.append( monotone_quintic_spline(*(sub_points.T), mids=1, ends=0) )
        print("adding..", end=" ", flush=True)
        p.add_func(f"Fit {i+1}", funcs[-1], [min(sub_points[:,0]),
                                             max(sub_points[:,1])], 
                   dash='dash', opacity=.3)
        print("done.", flush=True)

    avg_dist = lambda x: sum(min(1,max(0,f(x))) for f in funcs) / len(funcs)
    p.add_func(f"Avg", avg_dist, [min(x), max(x)], color=p.color(0))

    f = monotone_quintic_spline(*cdf_points(x), mids=1, ends=0)
    p.add_func(f"All", f, [min(x), max(x)], color=p.color(1))

    from scipy.stats import norm
    p.add_func(f"Truth", norm.cdf, [min(x), max(x)], color=p.color(3))

    p.show(y_range=[-1,2])
    exit()

