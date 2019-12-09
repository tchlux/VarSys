import numpy as np
# Local files.
from plot import Plot
from monotone import monotone_cubic_spline, monotone_quintic_spline


f = lambda x: np.sin(x) + x
df = lambda x: np.cos(x) + 1
ddf = lambda x: -np.sin(x)

num_points = 2

lower = 0
upper = 5/2 * np.pi


x = np.linspace(lower, upper, num_points)
y = f(x)
dy = df(x)
ddy = ddf(x)


cubic_fit = monotone_cubic_spline(x,values=np.vstack((y,dy)).T)
quintic_fit = monotone_quintic_spline(x,values=np.vstack((y,dy,ddy)).T)

legend = dict(
    xanchor = "center",
    yanchor = "top",
    x = .625,
    y = .205,
    orientation = "h",
    bgcolor="white",
    bordercolor="grey",
    borderwidth=.5
)


p = Plot("","","",font_family="times", font_size=20)
p.add("Nodes", [lower,upper], [f(lower), f(upper)])
p.add_func("sin(x) + x", f, [lower, upper], dash=None)
# p.add_func("df", df, [lower, upper])
# p.add_func("ddf", ddf, [lower, upper])
p.add_func("cubic", cubic_fit, [lower,upper], dash="dash")
p.add_func("quintic", quintic_fit, [lower,upper], dash="dot")
p.add("Nodes", [lower,upper], [f(lower), f(upper)], color=p.color(0), show_in_legend=False)
fig = p.plot(file_name="cubic-quintic-sin.html", legend=legend,
             width=700, height=400)

# import plotly
# help(plotly)
# help(plotly.offline.iplot) #(fig, show_link=False)
# plotly.offline.plot(fig, image="svg",
#                     filename="cubic-quintic-sin-test.pdf", 
#                     image_height=400, image_width=700)


# help(plotly.offline.plot)

# fig.write_image("cubic-quintic-test.pdf")

# offline.iplot({'data': [{'y': [4, 2, 3, 4]}], 
#                'layout': {'title': 'Test Plot'}},
#              image='svg')
