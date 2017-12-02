# Global variables
import sys
sys.path = ['/home/thomas/Git_Analytics/Analytics-Vis/Code'] + sys.path
sys.path = ['/home/thomas/Dropbox/Research/TL ML/Apriori_Tree'] + sys.path

import numpy as np

from settings import read_data_settings, data_for_spec, spec_str_to_dict, dict_to_spec_str
from tl_apriori_tree import AprioriTree


# Get the three primary command line arguments
group_id = "VirginiaTech"
project_name = ["VarSys", "VarSys_Per_Thread"]

# Bottom 20% of Variance
bot_20_var = "2.5e+08_&le;_V_&le;_7e+11"
# Mid 60% of Variance (20th - 80th percentile)
mid_60_var = "7e+11_&le;_V_&le;_7.7e+14"
# Top 20% of Variance
top_20_var = "7.7e+14_&le;_V_&le;_2.4e+17"

# # Top 10% of Variance
# top_10_var = "3.5e+15_&le;_V_&le;_2.4e+17"

var_settings = [(0.0, 20.0), (20.0, 80.0), (80.0, 100.0)]
var_setting_names = ("Bottom 20%", "Mid 60%", "Top 20%")

default_data_settings = read_data_settings(group_id, project_name[1])
data_settings_dict = dict(default_data_settings)
data_settings_order = [s[0] for s in default_data_settings]
# Update the settings that we want to lock down
data_settings_dict.update({
    "Mode":["Fread"],
    "Stdev":[],
    "Mean":[],
})

data_req = dict_to_spec_str(data_settings_dict, data_settings_order)

# Read the data from file to pass to plotting functions,
# change stdout so that the print statements from
# "data_for_spec" do not make their way to the terminal
normal_out = sys.stdout
sys.stdout = open("/dev/null",'w')
raw_data, header, types,_ = data_for_spec(
    group_id, project_name[1], spec_str_to_dict(data_req))
sys.stdout = normal_out

print("Raw data shape: (%s, %s)"%(raw_data.shape[0], len(header)))
print("Variance = [%s to %s]"%(min(raw_data["Variance"]), max(raw_data["Variance"])))
print()
print()
# Remove "Variance" and "Mode" attributes from header
header = [h for h in header if h not in ["Variance", "Mode"]]

perc_space = 8
print_col_widths = {h:2+max(len(str(v)) for v in [h]+data_settings_dict[h])
                    for h in header}
# Function for printing itemsets in a readable fashion
def pretty_itemset(iset):
    string = ""    
    values = {h:"--" for h in header}
    for i in iset:
        h, v = i.split(" -- ")
        values[h] = v
    for h in header:
        string += values[h] + " "*(print_col_widths[h] - len(values[h]))
    return string

print("Getting the natural occurrence of all length-1 items.")
data = raw_data[ header ]
data = [[h+" -- "+str(v) for h,v in zip(header,row)] for row in data]
min_support = 0.05
tree = AprioriTree(support=min_support)
original_item_sets = tree.mine_items(data, longest_item=1)
# Print out the itemsets
for size,item_set in enumerate(original_item_sets): 
    items = list(item_set.items())
    items.sort(key=lambda i: i[0])
    items.sort(key=lambda i: -i[1])

    print("  Item sets of length %i"%(size+1))
    print("="*25,"\n")

    print(" "*perc_space, ("%s"*len(header))%
          tuple(h+" "*(print_col_widths[h]-len(h)) for h in header))

    for i in items:
        perc_string = "%0.2f%%:"%(100.0*i[1])
        print(perc_string + " "*(perc_space - len(perc_string)),
              pretty_itemset(i[0]))
    print()
    print()

# Cycle through and generate itemset counts
for (lower, upper), setting_name in zip(var_settings, var_setting_names):
    # Calculate the actual numerical limits on variance given percentiles
    lower = np.percentile(raw_data["Variance"], lower)
    upper = np.percentile(raw_data["Variance"], upper)    
    # Reduce the raw data to just those rows with desird variance
    to_keep = [lower <= val <= upper for val in raw_data["Variance"]]
    data = raw_data[ np.where( to_keep ) ]

    # Print header for this setting
    print("="*(len(setting_name)+40))
    print(" "*20,setting_name)
    print("="*(len(setting_name)+40))
    print("Variance = [%s to %s]"%(lower,upper))
    print("Data shape: (%s,%s)"%(len(data), len(data.dtype.names)))
    print("Sample of data:")
    print(data[[np.random.randint(0,len(data)) for i in range(3)]])
    print()
    
    # Remove "Variance" and "Mode" attributes
    data = data[ header ]
    data = [[h+" -- "+str(v) for h,v in zip(header,row)] for row in data]

    # Assumed that data is a [[row1]...] where each row is unordered and
    # contains python types with an equals operator

    print("Generating Apriori Tree\n")

    min_support = 0.05
    tree = AprioriTree(support=min_support)
    item_sets = tree.mine_items(data)
        
    # Print out the itemsets
    for size,item_set in enumerate(item_sets): 
        if (size == 0):
            for i in item_set:
                if i in original_item_sets[0]:
                    item_set[i] = item_set[i] - original_item_sets[0][i]

        items = list(item_set.items())
        items.sort(key=lambda i: i[0])
        items.sort(key=lambda i: -i[1])

        print("  Item sets of length %i"%(size+1))
        print("="*25)
        if (size == 0): print("First column is percent shift in occurrence of the item set for this subset.")
        else: print("First column is the item set occurrence rate in this subset.")
        print()
        print(" "*perc_space, ("%s"*len(header))%
              tuple(h+" "*(print_col_widths[h]-len(h)) for h in header))

        for i in items:
            perc_string = "%0.2f%%:"%(100.0*i[1])
            print(perc_string + " "*(perc_space - len(perc_string)),
                  pretty_itemset(i[0]))
        print()
        print()
    
