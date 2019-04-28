import os
from util.data import Data
from util.system import load

configs, unique_configs = load("unique_configs.pkl")
print(configs)
print(len(unique_configs))
for uc in sorted(unique_configs):
    print(uc)


data_dir = "/Users/thomaslux/Git/LargeFiles/"

for data_file_name in sorted(os.listdir(data_dir)):
    if (data_file_name[:-4] not in {"2018-1", "2018-3"}): continue
    print("-"*70)
    print(data_file_name)
    try:
        file_path = os.path.join(data_dir,data_file_name)
        with open(file_path) as f:
            cols = f.readline().strip().split(",")
            config_cols = [c for c in cols if any(n[1:] in c.lower() for n in configs)]
        print(config_cols)
        print("Loading data..", flush=True)
        d = Data.load(file_path, sep=",", flush=True)
        print("Extracting unique configurations..", flush=True)
        d = d[config_cols]
        print("Getting unique values among those configs..", flush=True)
        d = d.unique()
        print("Mapping unique values into a set..", flush=True)
        contained_configs = set(map(tuple, d))
        print("Checking for content match..", flush=True)
        if any(c in contained_configs for c in unique_configs):
            print("This data may have the desired configurations..")
    except:
        print("Could not read..")
    print()


# print("Reducing to desired configs..")
# d = d[(i for i in range(len(d)) if
#        tuple(d[i,config_cols][0]) in unique_configs)]
