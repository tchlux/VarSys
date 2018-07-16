import numpy as np
from util.data import Data

from util.plot import Plot, multiplot


tm_file = "time_model_predicted[Delaunay].pkl"
print("Loading processed data..")
data = Data.load(tm_file)

print(data)

# min_pack = min(data[data.names[1]])
# data = data[data[data.names[1]] == min_pack]

def print_percentiles(data, percs=np.linspace(0,100,3), pretext=""):
    print(f"{pretext}{' ' if pretext else ''}Percentiles:")
    for p in percs:
        print(f"  {p:4.0f}th --", np.percentile(data,p))
    print(f" Range:", round(np.max(data) - np.min(data), 1))
    print()

first = True
for sys_config in data:
    l3, pack, time_model, pred_tm = sys_config
    # Make one plot
    name = f"Clock Cycles [L3 Size {l3} (KB), Packet Size {pack} (KB)]"
    p = Plot(name, "Byte Index", "Possible Value", "Clock Cycles")
    x = [v for v in time_model[time_model.names[0]]]
    y = [v for v in time_model[time_model.names[1]]]
    mean = time_model.names[2]
    mean_error = [t-p for (t,p) in zip(time_model[mean], pred_tm[mean])]
    name = f"[L3 {l3:04d}] [Packet {pack:04d}] Predicted mean absolute error"
    print_percentiles(list(map(abs,mean_error)), pretext=name)
    name = f"Predicted_Mean_Time_Model_[L3-{l3:04d}]_[Packet-{pack:04d}].csv"
    pred_tm.save(name)
    # # var = time_model.names[3]
    # p.add("Mean Error", x, y, mean_error, marker_line_width=.5, marker_line_color='#000')
    # p.plot(file_name="predicted_mean_var.html",
    #        append=(not first), show=(not first))
    first = False
