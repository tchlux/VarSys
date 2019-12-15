import numpy as np
# Include the "code" directory into the path and import code custom to this project.
import sys, os
sys.path += [os.path.abspath("code")]
from plot import Plot
sys.path.pop(-1)

c = "rgba(50,50,50,.5)"
a = .8
p = Plot("","α","β", font_family="times", font_size=18)
p.add("Region 1", [0, 3, 3], [3, 3, 0], mode="lines", fill="tonext",
      color=p.color(0), group=0)
p.add_func("Region 2", lambda x: (9 - x**2)**(1/2), [0,3], group=1,
           color=p.color(2), fill="tonext")
p.add("Region 3", [0, 3], [3, 0], mode="lines", fill="tonext", 
      color=p.color(3), group=2)
p.add("Region 4", [0, 1, 3], [3, 1, 0], mode="lines", fill="tozeroy", 
      color=p.color(1), group=3)

p.add("Point 1 to Region 1", [3,2], [4.5,3], color=c, group=0,
      line_color=p.color(0,alpha=a), mode="markers+lines", dash=None)
p.add("Point 2 to Region 2", [4,3*4/(4**2+3**2)**(1/2)], [3,3*3/(4**2+3**2)**(1/2)], 
      group=1, color=c, line_color=p.color(2, alpha=a),
      mode="markers+lines", dash="dot")
p.add("Point 3 to Region 3", [4,12/5], [1,3/5], color=c, group=2,
      line_color=p.color(3,alpha=a), mode="markers+lines", dash="dash")
p.add("Point 4 to Region 4", [1.5,9/13], [3.5,21/13], color=c, group=3,
      line_color=p.color(1,alpha=a), mode="markers+lines", dash="dashdot")

p.add_annotation(" 4 ", .3, .3, font_family="times", font_size=17,
                 show_arrow=False, border_width=1, x_anchor="left")
p.add_annotation(" 3 ", 1.02, 1.02, font_family="times", font_size=17,
                 show_arrow=False, border_width=1, x_anchor="left")
p.add_annotation(" 2 ", 1.55, 1.55, font_family="times", font_size=17,
                 show_arrow=False, border_width=1, x_anchor="left")
p.add_annotation(" 1 ", 2.3, 2.3, font_family="times", font_size=17,
                 show_arrow=False, border_width=1, x_anchor="left")


p.show(file_name="cubic_projection.html", width=400, height=400, show_legend=False)

