# ORDER:
#  [Frequency, File Size, Record Size, Threads, Iterations]

# Given a test that might be invalid, make it into the nearest valid test
def make_valid(test):
    # Record size, threads, and iterations should all be integers
    test[2] = round(int(test[2]))
    test[3] = round(int(test[3]))
    test[4] = round(int(test[4]))
    # Frequency should scale by the megahertz
    test[0] -= int(test[0])%100000
    # File size must be >= record size and a multiple of record size
    test[1] = max(test[2], test[1])
    test[1] = int(round(test[1] / test[2]) * test[2])

import numpy as np

lower = [2800000, 4, 4, 1, 200]
upper = [3500000, 16384, 8192, 64, 200]
steps = [8, 8, 8, 8, 1]
exp_spacing = [False, False, False, False, False]




values = []
for l,u,s,exp in zip(lower, upper, steps, exp_spacing):
    if exp:
        u = np.log2(u)
        l = np.log2(l)

    current_values = []
    if s > 1:
        for i in range(s):
            current_values += [ int(round(l + i*(u-l)/(s-1))) ]
    else:
        current_values += [(l + u) // 2]
        
    if exp:
        current_values = list(map(lambda v: 2**v, current_values))

    values.append( current_values )
        
# =============================
#      Hand Crafted Values     
# =============================

# Evenly spaced

# values = [
#     # Frequency
#     [2800000, 2900000, 3000000, 3100000, 3200000, 3300000, 3400000, 3500000],
#     # File Size
#     [4, 1024, 2048, 3072, 4096, 5120, 6144, 7168, 8192],
#     # Record Size
#     [4, 1024, 2048, 3072, 4096, 5120, 6144, 7168, 8192],
#     # Threads
#     [1, 10, 19, 28, 37, 46, 55, 64],
#     # Iterations
#     [200]
# ]

# Exponentially spaced

# values = [
#     # Frequency
#     [2800000, 2900000, 3000000, 3100000, 3200000, 3300000, 3400000, 3500000],
#     # File Size
#     [4, 16, 64, 256, 1024, 4096, ],
#     # Record Size
#     [4, 1024, 2048, 3072, 4096, 5120, 6144, 7168, 8192],
#     # Threads
#     [1, 10, 19, 28, 37, 46, 55, 64],
#     # Iterations
#     [200]
# ]

# values = np.meshgrid(*values)[0]
# print(values)
values = np.vstack(out.flatten() for out in np.meshgrid(*values)).T
# Make all of the tests valid (fix any problems with the full grid)
for test in values:
    make_valid(test)
# Now that tests are corrected, reduce to the set of unique tests
values = np.unique(values, axis=0)
# Sort them for consistency
sorted_indices = np.array(sorted(range(len(values)), key=lambda i: tuple(values[i])))
values = values[sorted_indices]

print()
names = ["Frequency", "File Size", "Record Size", "Threads", "Iterations"]
for n,v in zip(names,values.T):
    print(n,sorted(np.unique(v)))
print()

print("Number of tests:", len(values))

# Save them to a file that can be read by "estimate_runtime.py"
np.savetxt("grid.csv", values, fmt="%i", delimiter=",")

