import time
from scipy.spatial import Delaunay
import numpy as np
import pygal

def test(pts, dim):
    a = np.random.rand(pts, dim)*10
    start = time.time()
    hull = Delaunay(a)
    return time.time() - start

if __name__ == "__main__":
    # Should take about 5 minutes to run

    t_v_p = pygal.Line()
    time_vs_pts_2 = []
    time_vs_pts_3 = []
    time_vs_pts_4 = []
    time_vs_pts_5 = []
    now = time.time()
    for pts in range(0, 10001, 500):
        if pts == 0: pts = 100
        print("%.1f%% (%i points) at %.2f"%(
            100.0*pts/10001, pts,
            time.time() - now))
        now = time.time()
        time_vs_pts_2.append( test(pts, 2) )
        time_vs_pts_3.append( test(pts, 3) )
        time_vs_pts_4.append( test(pts, 4) )
        time_vs_pts_5.append( test(pts, 5) )
    t_v_p.add("Dim 5", time_vs_pts_5)
    t_v_p.add("Dim 4", time_vs_pts_4)
    t_v_p.add("Dim 3", time_vs_pts_3)
    t_v_p.add("Dim 2", time_vs_pts_2)
    t_v_p.title = "Delaunay Runtime"
    t_v_p.x_labels =  map(str,[100]+list(range(500, 10001, 500)))
    t_v_p.x_title = "Number of Points"
    t_v_p.y_title = "Computation Time (sec)"
    t_v_p.render_in_browser()
    t_v_p.render_to_file("Delaunay_time-vs-pts.svg")


    t_v_d = pygal.Line()
    time_vs_dim_20 = []
    time_vs_dim_50 = []
    time_vs_dim_80 = []
    time_vs_dim_100 = []
    now = time.time()
    dim = 2
    for i in range(9):
        print("(%i dimension) at %.2f"%(
            dim-1, time.time() - now))
        now = time.time()
        time_vs_dim_20.append(test(20, dim))
        time_vs_dim_50.append(test(50, dim))
        time_vs_dim_80.append(test(80, dim))
        time_vs_dim_100.append(test(100, dim))
        dim += 1
    print("(%i dimension) at %.2f"%(
        dim-1, time.time() - now))

    t_v_d.add("100 Points", time_vs_dim_100)
    t_v_d.add("80 Points", time_vs_dim_80)
    t_v_d.add("50 Points", time_vs_dim_50)
    t_v_d.add("20 Points", time_vs_dim_20)
    t_v_d.title = "Delaunay Runtime"
    t_v_d.x_labels = map(str, range(2,dim))
    t_v_d.x_title = "Number of Dimensions"
    t_v_d.y_title = "Computation Time (sec)"
    t_v_d.render_in_browser()
    t_v_d.render_to_file("Delaunay_time-vs-dim.svg")

