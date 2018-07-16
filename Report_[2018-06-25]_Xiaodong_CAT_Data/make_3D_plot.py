import numpy as np
from util.data import Data
from util.plot import Plot

_no_aesni = ""
# _no_aesni = "_no_aesni"
data_file_name = f"processed_study{_no_aesni}.dill"
print("Loading processed data..")
data = Data.load(data_file_name)
percentiles = np.linspace(0,1,11)
def percentile_color(p):
    low = np.array([50, 50,200])
    upp = np.array([200,50,50])
    return tuple(p*upp + (1-p)*low)

print(data[0,-1][0,-2])
print(data[0,-1][0,-1].inverse(percentiles))

min_pack = min(data[data.names[1]])
data = data[data[data.names[1]] == min_pack]

first = True
for sys_config in data:
    l3, pack, time_model = sys_config
    # Make one plot
    p = Plot(f"Distribution of Clock Cycles [L3 Size {l3} (KB), Packet Size {pack} (KB)]",
             "Byte Index", "Possible Value", "Clock Cycles")
    x = [v for v in time_model[time_model.names[0]]]
    y = [v for v in time_model[time_model.names[1]]]
    for i in percentiles:
        z = [cdf.inverse(i) for cdf in time_model[time_model.names[-1]]]
        color = p.color(color=percentile_color(i))
        p.add(f"{100*i:.0f}th Percentile", x, y, z, color=color,
              marker_line_width=.5, marker_line_color='#000')
        # for position in time_model:
        #     byte, val, clock_cycles, cdf = position
        #     print("",byte, val, clock_cycles)
    p.show(width=900, height=500, 
           file_name="clock_cycle_distributions.html" ,
           append=(not first), show=first)
    first = False
