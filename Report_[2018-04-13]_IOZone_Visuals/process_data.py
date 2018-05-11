import os
from util.data import read_struct, Struct
from util.decorators import cache


cur_dir = os.path.dirname(os.path.abspath(__file__))
files = sorted([os.path.join(cur_dir, f) 
                for f in os.listdir(cur_dir)
                if f[-4:] == ".csv"])

all_data_file = "data_full.pkl.gz"
small_data_file = "data_subset.pkl"
data = Struct()
use_full_data = True


if use_full_data:
    # Read the raw data if a pre-computed data file is not available.
    if not os.path.exists(all_data_file):
        types = (str, str, int, int, int, int, int, str, float, int)
        for path in files:
            print()
            print(f"Processing {path}..")
            data += read_struct(path, sep=",", types=types, verbose=True)
            print(data)
            print()
        data.save(all_data_file)
    else:
        # Otherwise, load the pre-computed data file.
        print("Loading full data...")
        data.load(all_data_file)

    print(data)
    # data.summarize(max_display=100)
    # sub_data = data[data["Test"] == "readers"]
    # print(sub_data)
    # sub_data.save("data_readers.pkl.gz")

    sub_data = data[data["Test"] == "random_readers"]
    sub_data.save("data_random_readers.pkl.gz")
    exit()
else:
    # Check for a small data file
    if not os.path.exists(small_data_file):
        # Read the raw data if a pre-computed data file is not available.
        if not os.path.exists(all_data_file):
            types = (str, str, int, int, int, int, int, str, float, int)
            for path in files:
                print()
                print(f"Processing {path}..")
                data += read_struct(path, sep=",", types=types, verbose=True)
                print(data)
                print()
            data.save(all_data_file)
        else:
            # Otherwise, load the pre-computed data file.
            print("Loading full data...")
            data.load(all_data_file)

        import random
        print("Generating random selection...")
        rows = list(range(len(data)))
        random.shuffle(rows)
        print("Extracting random rows...")
        small_data = data[rows[:10000]]
        print("Saving small data file...")
        small_data.save(small_data_file)
    else:
        print("Loading small data...")
        small_data = Struct().load(small_data_file)


    print(small_data)
    print()
    small_data.summarize()
