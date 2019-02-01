from util.data import Data
from util.stats import cdf_fit_func
from util.plot import Plot

d = Data.load("throughput_180913_Luna.csv")
d = d[:,:5] + Data(data=([row[5:]] for row in d), names=["Time (sec)"])
d.add_column(map(len,d["Time (sec)"]), name="Num Trials")
d.add_column(map(cdf_fit_func, d["Time (sec)"]), name="Time CDF")

print(d)

p = Plot(f"Time CDFs", "Time", "Probability a run has lower time")
for row in d:
    cache = row[d.names.index("Cache")]
    DB = row[d.names.index("DB")]
    p.add_func(f"{cache} Cache, {DB} DB", row[-1], row[-1]())
p.show(file_name="[2018-09-18]_time_CDFs_TF_data.html")
    
