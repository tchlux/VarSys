import os, random
from util.data import Data

# Set the random seed.
SEED = 0
random.seed(SEED)
SAMPLE_SIZE = 80

fname = "wendy-virtual.csv"
test = "readers"
data_name = fname + ".pkl"
sub_data_name = test + "_" + data_name


print("Loading data..")
if not os.path.exists(sub_data_name):
    print("Loading original data..")
    if os.path.exists(data_name):
        d = Data.load(data_name)
    else:
        d = Data.load(fname, sep=",", verbose=True)
        d.save(data_name)
        # 12002 (repeated header line)
        # 97498 (repeated header line)
        #  ... (4 more occurrences of this line)
    print("Done.")
    print()
    print(d)
    # Reduce to the "readers" test type
    d = d[d["Test"] == "readers"]
    print(d)
    # Collect by unique configurations
    params = ["Frequency", "File Size", "Record Size", "Num Threads"]
    d = d[params].unique().collect(d)
    d = d[params + ["Throughput"]]
    # Reduce all throughput sets with more than 80 samples to 80 samples.
    d["Throughput"] = (random.sample(tvals, SAMPLE_SIZE) for tvals in d["Throughput"])
    # Sort the data (to make it reasonably ordered)
    d.sort()
    print(d)
    # Save the data to file.
    d.save(sub_data_name)    
else:
    d = Data.load(sub_data_name)
print("Done.")
print()

# Create new data that mimics structure of old data.
nd = Data(
    names=["Machine","Store","Journal","Hyp","Hyp_Sched","VM_Sched","RCU","F_size","R_Size","Threads","Mode","Freq","Throughput"],
    types=[str,      str,    str,      str,  str,        str,       int,  int,     int,     int,      str,   int,   float]
)
for (freq, fs, rs, nt, thrpts) in d:
    # Cycle through each of the throughput values and add them to the new data.
    for thrpt in thrpts:
        nd.append(["new", "HDD", "yes", "xen", "CFQ", "NOOP", 128, fs, rs, nt, "Fread", freq, thrpt])

print(nd)
nd.save("MODIFIED-wendy-virtual.csv")
