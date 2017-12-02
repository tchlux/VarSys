import pygal

CUSTOM_CSS_FILE = 'file:///home/thomas/Git_Analytics/Analytics-Vis/Code/ViewData/Value_Distribution/no_dot.css'
num_points = 500
min_val = 0
width = 10

on_side = lambda x,knot,side: (x - knot) * side > 0
basis_1 = lambda x,knot,side=1: on_side(x,knot,side) * (x - knot)
basis_2 = lambda x,knot,side=1: on_side(x,knot,side) * (x - knot)**2
basis_3 = lambda x,knot,side=1: on_side(x,knot,side) * (x - knot)**3

funcs = [
    ("B-1, k-(4)",   lambda x: basis_1(x,4)),
    ("B-2, k-(4)",   lambda x: basis_2(x,4)),
    ("B-2, k-(4,2)", lambda x: basis_2(x,4) + basis_1(x,4,-1)),
    ("B-3, k-(2)",   lambda x: 0.01 * basis_3(x,2)),
]


# Create the plot
config = pygal.Config(fill=False, show_legend=True)
config.css.append(CUSTOM_CSS_FILE)
plot = pygal.XY(config)

for f_name, func in funcs:
    pts = []
    for i in range(num_points+1):
        x = (i / num_points) * width + min_val
        pts.append( (x, func(x)) )
    plot.add(f_name,pts)

plot.render_in_browser()
