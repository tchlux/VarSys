import numpy as np

#      data.py     
# =================
COMMON_SEPERATORS = [" ", ",", ";", "	"]

UPDATE_RATE = 10000  # How frequently to report progress percentages
MAX_ERROR_PRINTOUT = 10 # Only print out this many errors when processing data

NP_TYPES = {str:(np.str_, 16),    # Which numpy types correspond to 
            int:(np.int64,),      # the python types of variables
            float:(np.float64,)}
PY_TYPES = {value[0]:key for (key,value) in NP_TYPES.items()}

# Reduce a number to its minimal display form (no unnecessary 0's on right)
CLEAN_NUMBER_STRING = lambda number: str(number).rstrip("0").rstrip(".") \
                      if "." in str(number) else str(number)

# Print a numpy array space seperated so that it can be taken by fortran
def CLEAN_STRING(arr, string=""):
    for num in arr: string = string + " " + CLEAN_NUMBER_STRING(num)
    return string

CLEAN_ARRAY_STRING = lambda arr: (
    "\n".join([CLEAN_ARRAY_STRING(row) for row in arr]) if
    (len(arr.shape) > 1) else CLEAN_STRING(arr) )




# Return the type of a string (try converting to int and float)
def get_type(string):
    try:
        int(string)
        return int
    except ValueError:
        try:
            float(string)
            return float
        except ValueError:
            return str

# Pre:  "filename" is the name of an existing file (including path)
# Post: "filename" is scanned for viable seperators, those are all
#       characters used the same number of times on each line. If at
#       the end of the file there are still multiple candidate
#       seperators, they are sorted by most-occuring first. If one of
#       COMMON_SEPERATORS exist in the set, it is used, otherwise the
#       most occurring seperator is used.
def detect_seperator(filename, verbose=False):
    lines = []
    not_viable = set()
    with open(filename) as f:
        # Identify the potential seperators from the first line
        first_line = f.readline().strip()
        seperators = {char:first_line.count(char) for char in set(first_line)}
        # Cycle the file (narrowing down the list of potential seperators
        for line in f.readlines():
            line_chars = set(line.strip())
            # Make sure all tracked seperators appear the correct
            # number of times in the line
            for char in list(seperators.keys()):
                if (line.count(char) != seperators[char]):
                    not_viable.add( char )
                    seperators.pop( char )
            # WARNING: Do *not* break early, to verify that the
            #          seperators are used correctly throughout the
            #          entire file!
    # Get the seperators from the dictionary sorted by the decreasing
    # number of occurrences the character has on each line of the data.
    seperators = sorted(seperators.keys(), key=lambda k: -seperators[k])
    # Check to make sure there is a seperator detected
    if len(seperators) == 0:
        raise(Exception("ERROR: No consistently used seperator in file."))
    if verbose:
        print(" Detected the following possible seperators:\n  ", seperators)
    # Search through the common seperators first, using one of them if possible
    for sep in COMMON_SEPERATORS:
        if sep in seperators:
            if verbose: print(" Using '%s' as seperator."%sep)
            break
    else:
        # Otherwise, use the most occurring character that could be a seperator
        sep = max(seperators, key=lambda k: seperators[k])
        print("WARNING: Using '%s' as a seperator, even though it is unexpected."%(sep))
    # Return the automatically detected seperator
    return sep

# Pre:  "filename" is the name of a txt file that has a header
#       "sep" is the seperator used in the data file
#       "display" is True if the user wants the detected types printed
# Post: The header of the file as a list of strings,
#       the types of each header as a list,
#       the raw data as a [list of [lists of rows]]
def read_data_and_type(filename, sep=None, update_rate=UPDATE_RATE):
    if sep == None: sep = detect_seperator(filename, update_rate!=0)
    # Open the file to get the data
    with open(filename) as f:
        # Get the first line, assuming it's the header
        header = f.readline().strip().split(sep)
        # Get the first line of data and learn the types of each column
        line = f.readline().strip().split(sep)
        types = list(map(get_type, line))
        # Add the first line to the data store
        row_1 = tuple(t(v) for t,v in zip(types,line))
        # Get the rest of the raw data from the file and close it
        raw_data = f.readlines()
    # Print out a nicely formatted description of the detected types
    if update_rate != 0:
        print()
        print("Automatically detected types:")
        first_line  = ""
        second_line = ""
        for h,t in zip(header,types):
            first_line += "  "
            second_line += "  "
            length = max(len(h), len(t.__name__))
            first_line += h + " "*(length-len(h))
            second_line += t.__name__ + " "*(length-len(t.__name__))
        print(first_line)
        print(second_line)
    # Now process the rest of the data, dynamically overwriting old data
    to_remove = []
    for i,line in enumerate(raw_data):
        if (update_rate != 0) and not (i % update_rate):
            print("\r%0.1f%% complete"% (100.0*i/len(raw_data)), 
                  flush=True, end="")
        try:
            raw_data[i] = tuple(cast(value) for cast,value in
                                zip(types, line.strip().split(sep)))
        except:
            to_remove.append(i)
            if update_rate != 0:
                if len(to_remove) > MAX_ERROR_PRINTOUT:
                    if len(to_remove) == MAX_ERROR_PRINTOUT+1:
                        print("Suppressing output for errors... See all "+
                              "erroneous lines in seperate printout.")
                else:
                    print("\nERROR ON LINE %i, REMOVING FROM DATA:\n  %s\n"%(i+3,line))
    if update_rate != 0: print()
    # Remove lines that caused trouble when casting types
    for i in to_remove[::-1]:
        raw_data.pop(i)
    if (update_rate != 0) and (len(to_remove) > MAX_ERROR_PRINTOUT):
        print("ERRONEOUS LINES FROM DATA:\n  %s"%to_remove)
    # Return the header, types, and raw data
    return header, types, [row_1] + raw_data

# Given a file name, automatically detect the seperator and the types
# of each column, return all of this in a numpy structured array.
def read_struct(file_name, sep=None, verbose=False):
    print_rate = 0 if (not verbose) else UPDATE_RATE
    if verbose:
        print("Opening data file and reading into Python List...")
    try:
        header, types, data = read_data_and_type(file_name, sep, print_rate)
        # Sort the data by column precedent (standard python sort)
        # data.sort()
    except TypeError:
        msg = "Current version cannot read data with missing values\n"+\
              "            or values that change type throughout column. Types\n"+\
              "            are automatically detected based on first data row."
        raise(Exception("DATA ERROR: "+msg))
    # Convert the data into numpy form
    if verbose: print("Converting the data into numpy format...")
    dtypes = np.dtype([(h,)+NP_TYPES[t] for h,t in zip(header,types)])
    if verbose: print("Data types:",dtypes)
    data = np.array(data, dtype=dtypes)
    return data




# =======================
#      qHullDelaunay     
# =======================

# Wrapper class for using the Delaunay scipy code
class Delaunay:
    from scipy.spatial import Delaunay as Delaunay_via_qHull
    # For 1-D case, use regular interpolation
    from scipy import interpolate
    def __init__(self):
        self.pts = None
        self.values = None

    # Wrapper for 'predict' that returns a single value for a single
    # prediction, or an array of values for an array of predictions
    def __call__(self, x:np.ndarray, *args, **kwargs):
        single_response = len(x.shape) == 1
        if single_response:
            x = np.array([x])
        if len(x.shape) != 2:
            raise(Exception("ERROR: Bad input shape."))
        response = np.asarray(self.predict(x, *args, **kwargs), dtype=float)
        # Return the response values
        return response[0] if single_response else response

    # Use fortran code to compute the boxes for the given data
    def fit(self, x, y, normalize=True):
        self.pts = np.asarray(x.copy(), dtype=float)
        self.values = np.asarray(y.copy(), dtype=float)
        if normalize:
            # Normalize distances along each dimension
            self.shift = np.min(self.pts, axis=0)
            self.pts -= self.shift
            self.scale = np.max(self.pts, axis=0)
            self.pts /= self.scale
        else:
            self.shift = 0.0
            self.scale = 1.0
        # Handle special case for 1-D data, use regular linear
        #  interpolation, otherwise use standard Delaunay
        if (self.pts.shape[1] > 1):
            self.surf = Delaunay.Delaunay_via_qHull(self.pts)
        else:
            self.surf = Delaunay.interpolate.interp1d(self.pts[:,0],
                                                      self.values,
                                                      fill_value="extrapolate")

    # Use scipy code in order to evaluate delaunay triangulation
    def predict(self, x):
        # Compute the response values
        response = []
        for i,x_pt in enumerate(x):
            x_pt = (np.asarray(x_pt, dtype=float) - self.shift) / self.scale
            # Regular linear interpolation for 1-D data
            if (self.pts.shape[1] == 1):
                value = self.surf(x_pt[0])
                response.append(value)
                continue
            else:
                # Solve for the weights in the old Delaunay model
                simp_ind = self.surf.find_simplex(x_pt)
                # If a point is outside the convex hull, use the
                # closest simplex to extrapolate the value
                if simp_ind < 0: 
                    # Calculate the distance between the x_pt each simplex
                    simp_dists = self.surf.plane_distance(x_pt)
                    # Find the index of the closest simplex
                    simp_ind = np.argmin(abs(simp_dists))
                    
                # Solve for the response value
                simp = self.surf.simplices[simp_ind]
                #### Use the simplex from delauny to interpolate ####
                system = np.concatenate((self.pts[simp],
                    np.ones((simp.shape[0],1))), axis=1).T
                x_pt = np.concatenate((x_pt,[1]))
                try:
                    weights = np.linalg.solve(system, x_pt)
                except np.linalg.linalg.LinAlgError:
                    weights = np.linalg.lstsq(system[:-1,:], x_pt[:-1])[0]
                    weights /= np.sum(weights)
                    np.set_printoptions(formatter={"float":lambda f: "%.2f"%f})
                    # np.set_printoptions(precision=2)
                    print("ERROR")
                    print(system[:-1,:].T*self.scale + self.shift)
                    print(weights, np.sum(weights))
                    print(simp)
                    try: print("Distances:", simp_ind, simp_dists)
                    except: pass
                    print(x_pt)
                    exit()
                value = np.dot(self.values[simp],weights)
                #### Use Scipy's linear interpolator ###
                # plane = Delaunay.interpolate.LinearNDInterpolator(
                #     self.pts[simp], self.values[simp])
                # value = plane(x_pt)
                response.append(value)
        return response
