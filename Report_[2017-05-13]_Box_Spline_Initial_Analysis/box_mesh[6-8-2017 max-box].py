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
def compute_box(x, center=None, low_width=None, upp_width=None, func=quadratic):
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

    # Evaluates the underlying box function for this box at a given point
    def __call__(self, x):
        # Return the computed quadratic box spline evaluation
        return compute_box(x, self.center, self.low_width*self.width_scalar,
                           self.upp_width*self.width_scalar, self.box_func) 


# Function for creating a surface over a region defined by the
# axis-aligned bounding box of all data points
class BoxMesh:
    def __init__(self, box_func=quadratic, width_scalar=1.0):
        self.boxes = []
        self.values = []
        self.box_func = box_func
        self.width_scalar = width_scalar

    # Given a set of control points, this constructs the maximum box
    # mesh, that is to say the boxes around each control point have
    # the largest minimum side length possible.
    def add_points(self, control_points, max_box=True):
        self.values = control_points[:,-1]
        control_points = control_points[:,:-1]

        # Cycle through the box centers and construct the max-box
        # around each point using the dimension-sorted data
        for i,(center,value) in enumerate(zip(control_points,values)):
            box = Box(center, func=self.box_func, width_scalar=self.width_scalar)

            # Holder for candidate boundary points for this box
            candidates = list(range(len(control_points)))
            candidates.pop(i)
            candidates.sort(key = lambda c: np.max(abs(center - control_points[c])))
            
            # Holder for the neighbors of this box (on its edge)
            neighbors = []

            # Cycle through and shrink the box, touching immediate neighbors
            for c in candidates:
                if not box.contains(control_points[c]): continue

                other_pt = control_points[c]
                # handle choosing a split dimension normally
                split_dim = np.argmax(abs(center - other_pt))

                # Shrink 'box' along that best dimension
                width = center[split_dim] - other_pt[split_dim]
                if width < 0:
                    box.upp_width[split_dim] = min(box.upp_width[split_dim],-width)
                else:
                    box.low_width[split_dim] = min(box.low_width[split_dim], width)
                
                # Update the neighbors for the new box
                neighbors.append( [split_dim,c] )

            # Skip the rest of the computations if we don't need max boxes
            if not max_box:
                self.boxes.append(box)
                continue

            #     Constructing the Box with Largest Minimum Side Length     
            #===============================================================

            # Temporary variable names
            dim = control_points.shape[1]
            num_pts = control_points.shape[0]
            min_bounds = np.min(control_points, axis=0)
            max_bounds = np.max(control_points, axis=0)

            # Try and find a larger minimum side length, once no
            # improvement can be made, we must have the max-box
            for _ in range(num_pts):
                # Adjust the side lengths, ignoring infinities
                low_width = np.where(box.low_width < float('inf'),
                                     box.low_width, box.center - min_bounds)
                upp_width = np.where(box.upp_width < float('inf'),
                                     box.upp_width, max_bounds - box.center)
                side_lengths = low_width + upp_width
                # Get the smallest side and its length
                min_side = np.argmin(side_lengths)
                min_len = np.min(side_lengths)
                # Get the current neighbors on this minimum side
                to_use = [n for n in neighbors if n[0] == min_side]

                for d in range(dim):
                    # Don't try and shorten the min side
                    if d == min_side: continue
                    # Holder for knowing if we've improved the box
                    improved = False
                    
                    # At most 2 on non-gridded data
                    for pt in to_use:
                        # Get this neighboring point
                        neighbor_pt = control_points[pt[1]]

                        # Shorten the existing box on dimension d to stop at this neighbor
                        if neighbor_pt[d] < center[d]:
                            old_width = (box.low_width, box.low_width[d])
                            box.low_width[d] = center[d] - neighbor_pt[d]
                            old_neighbor = [n for n in neighbors if (
                                n[0] == d and control_points[n[1],d] < center[d])]
                        else:
                            old_width = (box.upp_width, box.upp_width[d])
                            box.upp_width[d] = neighbor_pt[d] - center[d]
                            old_neighbor = [n for n in neighbors if (
                                n[0] == d and control_points[n[1],d] > center[d])]

                        # Extend the box on min_side to infinity
                        if neighbor_pt[min_side] < center[min_side]:
                            growing = box.low_width
                        else:
                            growing = box.upp_width

                        old_min_side = growing[min_side]
                        growing[min_side] = float('inf')
                        new_neighbor = None

                        new_neighbor = None
                        # Identify the new boundary for the box
                        for c in candidates:
                            if box.contains(control_points[c]):
                                if ((new_neighbor == None) or
                                    ( abs(center[min_side] - control_points[c][min_side])
                                      <
                                      abs(center[min_side] - control_points[new_neighbor[1]][min_side])
                                    )):
                                    growing[min_side] = abs(center[min_side] - control_points[c][min_side])
                                    new_neighbor = [min_side,c]

                        # Adjust the side lengths, ignoring infinities
                        low_width = np.where(box.low_width < float('inf'),
                                             box.low_width, box.center - min_bounds)
                        upp_width = np.where(box.upp_width < float('inf'),
                                             box.upp_width, max_bounds - box.center)
                        side_lengths = low_width + upp_width

                        # Identify the new minimum dimension of the box
                        if np.min(side_lengths) > min_len:
                            # If we've increased the minimum side length
                            improved = True
                            # Formally adjust the modified neighbor boundary
                            pt[0] = d
                            # Add the new neighbor that was encountered
                            if new_neighbor:
                                neighbors.append(new_neighbor)
                            # Remove the old neighbor on this dimension
                            for to_remove in old_neighbor:
                                neighbors.remove(to_remove)
                            break
                        else:
                            # Otherwise, revert changes
                            old_width[0][d] = old_width[1]
                            growing[min_side] = old_min_side
                            pt[0] = min_side

                    # Improvement found, break out of this iteration
                    if improved: break

                else:
                    # If we've checked all dimensions and no
                    # improvement could be found, then we are done!
                    break

            self.boxes.append(box)

        # Find the value at the center of a box (same for all boxes)
        self.box_center_value = self.boxes[0](self.boxes[0].center)

    # Evaluate all of the box splines and return the final function
    # value at the requested point.
    def __call__(self, x):
        # If we're using a linear function, calculate lagrange-style interpolation value
        if self.box_func == linear:
            value = 0
            for (b,v) in zip(self.boxes, self.values):
                box_val = b(x)
                value += box_val * v / self.box_center_value
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
            return None

        # Return the evaluation of all box splines for this point
        return numerator / denominator


# =======================
#      Plotting Code     
# =======================

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
    plot.add("%s boundary"%(box.id), *list(zip(*corners)),
             mode='lines', color=color, opacity=opacity)
    return color


# Draw the boundary boxes for a given mesh
def draw_boxes_2D(plot, mesh):
    extra = 0.1
    min_max_x = (min(b.center[0] for b in mesh.boxes)-extra, max(b.center[0] for b in mesh.boxes)+extra)
    min_max_y = (min(b.center[1] for b in mesh.boxes)-extra, max(b.center[1] for b in mesh.boxes)+extra)
    min_max = list(zip(min_max_x, min_max_y))
    colors = []
    # First, draw the borders of all of the boxes
    for box in mesh.boxes:
        colors.append( draw_2D_box(plot, box, min_max) )

    # Draw the centers second in order to make them on 'top' of the borders
    for box,color in zip(mesh.boxes, colors):
        plot.add("%s center"%(box.id), *[[v] for v in box.center],
                 show_in_legend=False, marker_size=5, symbol='square',
                 color=color)


# ====================================================
#      Testing the code for generating box meshes     
# ====================================================

if __name__ == "__main__":

    use_max_box = False
    random_points = False
    num_points = 8

    func = linear
    size = 1.0
    plot_range = [-0.2, 1.2]
    plot_points = 3000

    # fun = lambda x: (x[0]-num_points/3)*(x[1]-num_points/2)**2
    fun = lambda x: np.cos(x[0]) * np.sin(x[1])

    # Generate testing points in a grid pattern
    points = np.meshgrid(range(num_points), range(num_points))
    points = np.array([points[0].flatten(), points[1].flatten()]).T

    if random_points:
        # Generate testing points randomly spaced
        points = np.random.random(size=(num_points**2,2)) * num_points

    # Pre-defined random points for re-use in testing
    # points = np.array([
    #     [ 1.12392062,  2.8193117 ],
    #     [ 2.41967177,  1.85301683],
    #     [ 2.13690637,  0.61453008],
    #     [ 1.83613747,  2.78149679],
    #     [ 1.44421677,  0.31735115],
    #     [ 0.19415758,  1.48313631],
    #     [ 2.96359513,  2.83153523],
    #     [ 2.84295907,  0.99043686],
    #     [ 0.36184628,  1.51468645],
    #     [ 2.69452803,  1.58899823],
    #     [ 2.05189331,  1.93667706],
    #     [ 0.08524652,  2.59040647],
    #     [ 2.63828929,  0.71555996],
    #     [ 0.96454286,  2.08487523],
    #     [ 0.6620956,   2.81263217],
    #     [ 1.61650631,  0.58805694],
    # ])

    # Calculate the associated response values
    values = np.array([[fun(pt) for pt in points]]).T
    points = np.concatenate((points, values), axis=1)
    points = points[points[:,1].argsort()]
    points = points[points[:,0].argsort()]

    # print(points)

    # Normalize the points
    points[:,2] -= min(points[:,-1])
    max_val = max(max(r[:-1]) for r in points)
    points[:,:2] /= max_val
    points[:,2] += 1.3

    print("Constructing mesh...")
    surf = BoxMesh(func, size)
    surf.add_points(points, use_max_box)
    print("Creating visualization...")
    surf_p = Plot()
    print(" adding control points...")
    surf_p.add("Control Points", *(points.T))
    print(" adding mesh surface...")
    surf_p.add_func("Surface", surf, plot_range, plot_range,
                    use_gradient=True, plot_points=plot_points)

    surf_p.plot()

    # print("Drawing plots...")
    # boxes_p = Plot()
    # draw_boxes_2D(boxes_p, surf)
    # multiplot([[boxes_p,surf_p]])

    # print("Constructing second mesh...")
    # surf2 = BoxMesh(func, size)
    # surf2.add_points(points, not use_max_box)
    # print("Drawing plots...")
    # boxes_p2 = Plot()
    # draw_boxes_2D(boxes_p2, surf2)
    # surf_p2 = Plot()
    # surf_p2.add("Control Points", *(points.T))
    # surf_p2.add_func("Surface", surf2, plot_range, plot_range, use_gradient=True)


    # multiplot([[boxes_p,surf_p],
    #            [boxes_p2,surf_p2]])


