import numpy as np

# Pre:  "x" is a 1 x N dimensional numpy array or list in the range [0,1]
# Post: The value of the quadratic box spline in N dimensions
#       is computed for the given point.
def linear(x):
    x *= 2
    # Compute the box spline function value for the x-coordinate based
    # on which region of the box that it is in
    func_val = 1.0
    for d in range(len(x)):
        if (0 <= x[d] < 1):
            func_val *= x[d]
        elif (1 <= x[d] < 2):
            func_val *= (2 - x[d])
    return func_val

# Pre:  "x" is a 1 x N dimensional numpy array or list in the range [0,1]
# Post: The value of the quadratic box spline in N dimensions
#       is computed for the given point.
def quadratic(x):
    x *= 3
    # Compute the box spline function value for the x-coordinate based
    # on which region of the box that it is in
    func_val = 1 / 2**(x.shape[0])
    for d in range(len(x)):
        if (0 <= x[d] < 1):
            func_val *= x[d]**2
        elif (1 <= x[d] < 2):
            func_val *= -(2*x[d]**2 - 6*x[d] + 3)
        elif (2 <= x[d] <= 3):
            func_val *= (x[d] - 3)**2
    return func_val

# Pre:  "x" is a 1 x N dimensional numpy array or list in the range [0,1]
# Post: The value of the classic cubic box spline in N dimensions
#       is computed for the given point.
def cubic(x):
    x *= 4
    # Compute the box spline function value for the x-coordinate based
    # on which region of the box that it is in
    func_val = 1 / (2*3**(x.shape[0]))
    for d in range(len(x)):
        if (0 <= x[d] < 1):
            func_val *= (x[d]**3)
        elif (1 <= x[d] < 2):
            func_val *= (4 - 12*x[d] + 12*x[d]**2 - 3*x[d]**3)
        elif (2 <= x[d] < 3):
            func_val *= (-44 + 60*x[d] - 24*x[d]**2 + 3*x[d]**3)
        elif (3 <= x[d] <= 4):
            func_val *= (4-x[d])**3
    return func_val

# Uses a guassian distribution that doesn't go to zero, but rather 3
# standard deviations at boundaries
def gaussian(x):
    x -= 0.5
    x *= 7
    func_val = np.prod(np.exp(-x**2/2) / np.sqrt(2*np.pi))
    return func_val

# Pre:  "x" is a 1 x N dimensional numpy array or list
#       "center" is a 1 x N dimensional numpy array or list
#       representing the bottom left (in cartesian coordinates) of the box
#       "[.*]_width" is the set of widths of this box spline in the
#       less than center and greater than center directions
#       "func" is the box function to use assuming we have x scaled
#       into the unit cube
# Post: The value of the classic quadratic box spline in N dimensions
#       is computed for the given point.
def compute_box(x, center=None, low_width=None, upp_width=None, func=linear):
    # Initialize center and width of box 3
    if type(center) == type(None):
        center = np.ones(x.shape)/2
    if type(low_width) == type(None):
        low_width = np.ones(x.shape)/2
    if type(upp_width) == type(None):
        upp_width = np.ones(x.shape)/2

    # Make sure we don't modify the input x value when evaluating
    x = x.copy()
    # Scale the point to be in the space where the box is the unit cube
    x -= center
    x = np.where(x < 0, (1 - x / low_width)/2, (1 + x / upp_width)/2)
    # If the point is outside of the box, return 0
    if (min(x) <= 0) or (max(x) >= 1): return 0

    # Return the final function value at that point
    return func(x)
    

# Function for storing a center and width associated with a quadratic
# box spline as well as evaluating the box spline
class Box:
    box_count = -1 # Used for generating ID numbers
    # Organize the necessary values for the box
    def __init__(self, center, low_width=None, upp_width=None,
                 func=linear, width_scalar=1.0):
        # Initialize lower and upper widths if they are not given
        if type(low_width) == type(None):
            low_width = np.ones(center.shape) * float('inf')
        if type(upp_width) == type(None):
            upp_width = np.ones(center.shape) * float('inf')
        # Get an id for this box
        self.id = Box.box_count = Box.box_count + 1
        self.box_func = func
        self.width_scalar = width_scalar
        self.center = center
        self.low_width = low_width
        self.upp_width = upp_width

    # Returns true if this box contains "pt"
    def contains(self, pt):
        return (np.all(self.center-self.low_width < pt) and
                np.all(self.center+self.upp_width > pt))

    # String representation of this box
    def __str__(self):
        string =  "ID: %s\n"%self.id
        string += "  center: %s\n"%self.center
        string += "  lower:  %s\n"%(self.center - self.low_width)
        string += "  upper:  %s\n"%(self.center + self.upp_width)
        return string

    # Evaluates the underlying box function for this box at a given point
    def __call__(self, x):
        # Return the computed quadratic box spline evaluation
        return compute_box(x, self.center, self.low_width*self.width_scalar,
                           self.upp_width*self.width_scalar, self.box_func) 


# Function for creating a surface over a region defined by the
# axis-aligned bounding box of all data points
class MaxBoxMesh:
    def __init__(self, box_func=linear, width_scalar=1.0):
        self.boxes = []
        self.values = []
        self.box_func = box_func
        self.width_scalar = width_scalar

    # Given a set of control points, this constructs the maximum box
    # mesh, that is to say the boxes around each control point have
    # the largest minimum side length possible.
    # Computational complexity: O(n^2*log(n) * d)
    def fit(self, control_points, values):
        self.boxes = []
        self.values = values

        # Cycle through the box centers and construct the max-boxes
        for i,(center,value) in enumerate(zip(control_points,values)):
            box = Box(center, func=self.box_func, width_scalar=self.width_scalar)

            # Holder for candidate boundary points for this box
            candidates = list(range(len(control_points)))
            candidates.pop(i)
            candidates.sort(key = lambda c: np.max(abs(center - control_points[c])))

            # Cycle through and shrink the box, touching immediate neighbors
            for c in candidates:
                if not box.contains(control_points[c]): continue

                other_pt = control_points[c]
                # handle choosing a split dimension (max distance axis)
                split_dim = np.argmax(abs(center - other_pt))

                # Shrink 'box' along that best dimension
                width = center[split_dim] - other_pt[split_dim]
                if width < 0:
                    box.upp_width[split_dim] = min(box.upp_width[split_dim],-width)
                else:
                    box.low_width[split_dim] = min(box.low_width[split_dim], width)

            # Now we have found the box with the largest minimum
            # distance between the center and its boundary
            self.boxes.append(box)

        # Store the value at the center of a box (same for all boxes)
        self.box_center_value = self.boxes[0](self.boxes[0].center)
        # Calculate the lipschitz constant for this data
        self.calculate_lipschitz()

    # Stores the lipschitz constant of this data set in self.lipschitz
    # Computational complexity: O(n^2)
    def calculate_lipschitz(self):
        self.lipschitz = 0
        # Calculate the lipschitz constant of the data
        for (b,v) in zip(self.boxes,self.values):
            for (other_b,other_v) in zip(self.boxes,self.values):
                if (v == other_v): continue
                dist = np.sqrt(np.sum((b.center - other_b.center)**2))
                self.lipschitz = max(self.lipschitz, abs(v - other_v) / dist)
        return self.lipschitz

    # Calculate and return the plus/minus lipschitz error estimate at a given x-value
    def error(self, x):
        # Get the surface value
        surf_val = self(x)
        # Initialize holder for plus and minus error
        error = [float('inf'),float('inf')]
        for (b,v) in zip(self.boxes, self.values):
            dist = np.sqrt(np.sum( (b.center - x)**2 ))
            error[0] = min(error[0], abs(v - dist*self.lipschitz - surf_val))
            error[1] = min(error[1], abs(v + dist*self.lipschitz - surf_val))
        # Remove all infinities from error
        while float('inf') in error:
            error[error.index(float('inf'))] = 0
        # Return the plus and minus error at a given data point
        return error

    # Evaluate all of the box splines and return the final function
    # value at the requested point.
    def __call__(self, x):
        # Make sure x is in an array
        if len(x.shape) == 1:
            x = np.array([x])
        if len(x.shape) != 2:
            raise(Exception("ERROR: Bad input x shape."))
        # Use the NURBS approach to smoothing
        numerator = np.zeros(x.shape[0])
        denominator = np.zeros(x.shape[0])
        # Calculate the numerator and denominator
        for (b,v) in zip(self.boxes, self.values):
            box_vals = np.array(tuple(b(pt) for pt in x))
            numerator += box_vals * v
            denominator += box_vals

        # Adjust for points outside of the region
        if 0 in denominator:
            print("WARNING: Some points are out of the range of all boxes.")

        # Return the evaluation of all box splines for this point
        response = np.where(denominator != 0, numerator / denominator, float('inf'))
        return response[0] if len(response) == 1 else response
