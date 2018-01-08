import numpy as np
from plotly_interface import *
import os

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
    box_count = 0

    def __init__(self, center, low_width=None, upp_width=None,
                 func=quadratic, width_scalar=1.0):
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

    def __call__(self, x):
        # Return the computed quadratic box spline evaluation
        return compute_box(x, self.center, self.low_width*self.width_scalar,
                           self.upp_width*self.width_scalar, self.box_func) 


# Function for creating a surface over the region defined by the
# convex hull of a set of control points.
class BoxMesh:
    def __init__(self, box_func=quadratic, width=1.0):
        self.boxes = []
        self.values = []
        self.box_func = box_func
        self.width_scalar = width

    # Function for changing the widths of internally kept boxes based
    # on the location of the added control point, maximizing minimum side length
    def add_point(self, control_point):
        # Extract the value from the control point and create a new box
        value = control_point[-1]
        point = control_point[:-1]
        # Identify the existing boxes that contain this new point
        boxes = [b for b in self.boxes if b.contains(point)]
        # If there are no boxes, then create a new one
        if len(boxes) == 0:
            widths = np.ones(shape=point.shape) * float('inf')
            new_box = Box(point, widths.copy(), widths.copy(),
                          self.box_func, self.width_scalar,
                          len(self.boxes))
        else:
            # Generate a new box with the largest possible dimensions,
            # knowing that this will cause it to overlap with all the
            # existing boxes that contain it
            low_width =  point - np.min([(box.center - box.low_width) for box in boxes], axis=0)
            upp_width = -point + np.max([(box.center + box.upp_width) for box in boxes], axis=0)
            new_box = Box(point, low_width, upp_width, self.box_func, self.width_scalar)

            # Cycle through and shrink the box however is necessary in
            # order to make it not contain the other boxes while
            # keeping its size as large as possible
            for box in boxes:
                # Identify dimension with the largest gap
                split_dim = np.argmax(abs(box.center - new_box.center))
                # Split along that dimension
                width = abs(box.center[split_dim] - new_box.center[split_dim])
                if box.center[split_dim] < new_box.center[split_dim]:
                    lower, upper = box, new_box
                else:
                    lower, upper = new_box, box
                # Update the boundaries between the two boxes
                lower.upp_width[split_dim] = min(lower.upp_width[split_dim],width)
                upper.low_width[split_dim] = min(upper.low_width[split_dim],width)
            
        # Store the control point value and the new box
        self.values.append(value)
        self.boxes.append(new_box)

    # Evaluate all of the box splines and return the final function
    # value at the requested point.
    def __call__(self, x):
        # If we're using a linear function, calculate lagrange-style interpolation value
        if self.box_func == linear:
            value = 0
            for (b,v) in zip(self.boxes, self.values):
                box_val = b(x)
                value += box_val * v / b(b.center)
            return value

        # Otherwise, use the NURBS approach to smoothing
        numerator = denominator = 0
        # Calculate the numerator and denominator
        for (b,v) in zip(self.boxes, self.values):
            box_val = b(x)
            numerator += box_val * v
            denominator += box_val

        # Adjust for points outside of the region
        if denominator == 0:
            print("WARNING: Out of the range of all boxes.")
            numerator = 0
            denominator = 1

        # Return the evaluation of all box splines for this point
        return numerator / denominator

# Produces a list of the points that define the corners of the given
# box, when the width is infinity, uses min_max
def box_corners(box, low_upp):
    points = []
    bits = [0] * len(box.center)
    for d in range(2*len(bits)):
        bit_index = d % len(bits)
        bits[bit_index] = (bits[bit_index] + 1) % 2
        # Scale infinite sides to be just wide of the boundary points
        low_width = np.where(abs(box.low_width) != float('inf'), box.low_width, low_upp[0])
        upp_width = np.where(abs(box.upp_width) != float('inf'), box.upp_width, low_upp[1])
        # Add the next corner
        points += [[ box.center[i] + (-low_width[i] if bits[i] else upp_width[i])
                     for i in range(len(bits)) ]]
    return points

# Draw a box in a 2D plotly plot
def draw_2D_box(plot, box, min_max):
    # Generate the absolute lower and upper plotting bounds
    low_upp = [[box.center[i] - min_max[0][i],
                min_max[1][i] - box.center[i]] for i in range(len(box.center))]
    low_upp = np.array(low_upp).T
    # Get the boundary points of the box
    corners = box_corners(box, low_upp)
    corners.append(corners[0])
    corners = np.array(corners)
    # Add the box to the plot
    opacity = 0.7
    plot.color_num += 1
    color = plot.color(plot.color_num, alpha=opacity)
    center = list(box.center)
    plot.add("%s center"%(box.id), *[[v] for v in center],
             group=str(box.id), marker_size=5, symbol='square', color=color)

    try:
        plot.add("%s boundary"%(box.id), *list(zip(*corners)),
             group=str(box.id), show_in_legend=False,
             mode='lines', color=color, opacity=opacity)
    except:
        print(corners)


# Draw the boundary boxes for a given mesh
def draw_boxes_2D(plot, mesh):
    extra = 0.1
    min_max_x = (min(b.center[0] for b in mesh.boxes)-extra, max(b.center[0] for b in mesh.boxes)+extra)
    min_max_y = (min(b.center[1] for b in mesh.boxes)-extra, max(b.center[1] for b in mesh.boxes)+extra)
    min_max = list(zip(min_max_x, min_max_y))
    for box in mesh.boxes:
        draw_2D_box(plot, box, min_max)


# Create a minimal box mesh over a set of points given an error tolerance
def minimal_mesh(points, func=linear, size=1.0, rel_error_tol=0.01, show_steps=False):
    surf = BoxMesh(func, size)

    # Start the mesh off with the median valued point
    vals = [(v,i) for (i,v) in enumerate(points[:,-1])]
    vals.sort()
    median_point = points[ vals[len(vals)//2][1] ]
    points = np.delete(points, vals[len(vals)//2][1], 0 )
    surf.add_point(median_point)

    if show_steps:
        os.system("echo '' > steps.html")
        p = Plot()
        draw_boxes_2D(p, surf)
        p.plot(show=False)
        os.system("cat temp-plot.html >> steps.html")


    # Calculate the initial relative error.
    rel_error = abs(np.array([surf(pt[:-1]) for pt in points]) - points[:,-1]) / np.where(points[:,-1]==0,0.1,points[:,-1])
    # Add points until the maximum relative error is below some threshold
    for i in range(len(points)):
        if max(rel_error) < rel_error_tol: break
        pt_to_add = points[np.argmax(rel_error)]
        points = np.delete(points, np.argmax(rel_error), 0 )
        surf.add_point( pt_to_add )
        rel_error = abs(np.array([surf(pt[:-1]) for pt in points]) - points[:,-1]) / np.where(points[:,-1]==0,0.1,points[:,-1])

        if show_steps:
            p = Plot()
            draw_boxes_2D(p, surf)
            p.plot(show=False)
            os.system("cat temp-plot.html >> steps.html")

    if show_steps: os.system("open steps.html")

    return surf



# Testing the code for generating box spline meshes
if __name__ == "__main__":

    dist = 2
    fun = lambda x: (x[0]-dist/3)*(x[1]-dist/2)**2
    # fun = lambda x: np.cos(x[0]) * np.tan(x[1])

    # Generate testing points in a grid pattern
    points = np.meshgrid(range(dist+1), range(dist+1))
    points = np.array([points[0].flatten(), points[1].flatten()]).T
    values = np.array([[fun(pt) for pt in points]]).T
    points = np.concatenate((points, values), axis=1)
    points[:,:2] /= dist
    points[:,2] -= min(points[:,2])
    points[:,2] += 1.3

    print(points)


    minimal_mesh(points, show_steps=True)


    exit()

    # =======================================
    #      Simple example with two plots     
    # =======================================

    smoothness = 1.0

    # =========================================
    #      Linear Box Spline Interpolation     
    # =========================================

    p = Plot()
    p.add("True Points", points[:,0], points[:,1], points[:,2])
    plot_range = [0,1]

    print("Adding surface...")
    surf = minimal_mesh(points, func=linear, size=smoothness)

    p.add_func("Linear Surface", lambda x: surf(x), plot_range,
               plot_range, opacity=0.7, use_gradient=True) 

    # print()
    # for b in surf.boxes:
    #     p.add_func("Box "+str(b.id)+" weight", lambda x: b(x), plot_range,
    #                plot_range, mode='markers', marker_size=3,
    #                marker_line_width=1)
    #     print(b)

    print("Generating html...")
    p.plot("Linear Box Spline Interpolation",fixed=False,show=False)

    os.system("cat temp-plot.html > temp.html")

    # ============================================
    #      Quadratic Box Spline Interpolation     
    # ============================================

    p = Plot()
    p.add("True Points", points[:,0], points[:,1], points[:,2])
    plot_range = [0,1]

    print("Adding surface...")
    surf = minimal_mesh(points, func=quadratic, size=smoothness)
    p.add_func("Quadratic Surface", lambda x: surf(x), plot_range,
               plot_range, opacity=0.7, use_gradient=True) 

    # print()
    # for b in surf.boxes:
    #     p.add_func("Box "+str(b.id)+" weight", lambda x: b(x), plot_range,
    #                plot_range, mode='markers', marker_size=3,
    #                marker_line_width=1)
    #     print(b)

    print("Generating html...")
    p.plot("Quadratic Box Spline Interpolation",fixed=False,show=False)

    os.system("cat temp-plot.html >> temp.html")

    os.system("open temp.html")
