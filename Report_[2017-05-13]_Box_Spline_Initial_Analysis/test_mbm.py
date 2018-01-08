from max_box_mesh import *
from plotly_interface import *

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
        low_width = np.where(abs(box.low_width) != float('inf'), box.width_scalar*box.low_width, low_upp[0])
        upp_width = np.where(abs(box.upp_width) != float('inf'), box.width_scalar*box.upp_width, low_upp[1])
        # Add the next corner
        points += [[ box.center[i] + (-low_width[i] if bits[i] else upp_width[i])
                     for i in range(len(bits)) ]]
    return points

# Draw a box in a 2D plotly plot
def draw_box(plot, box, min_max):
    # Generate the absolute lower and upper plotting bounds
    low_upp = [[box.center[i] - min_max[0][i],
                min_max[1][i] - box.center[i]] for i in range(len(box.center))]
    low_upp = np.array(low_upp).T
    # Get the boundary points of the box
    corners = box_corners(box, low_upp)
    corners.append(corners[0])
    corners = np.array(corners)
    # If we're one dimensional, force y coordinate to be zero
    if corners.shape[1] == 1:
        corners = np.concatenate((corners, np.zeros(shape=corners.shape)), axis=1)
    # Add the box to the plot
    opacity = 0.7
    plot.color_num += 1
    color = plot.color(plot.color_num, alpha=opacity)
    center = list(box.center)
    plot.add("%s boundary"%(box.id), *list(zip(*corners)),
             mode='lines', color=color, opacity=opacity)
    return color


# Draw the boundary boxes for a given mesh
def draw_boxes(plot, mesh):
    centers = np.array([b.center for b in mesh.boxes])
    min_max = np.array([np.min(centers, axis=0), np.max(centers, axis=0)])
    extra = 0.1 * (min_max[1,:] - min_max[0,:])
    min_max[0,:] -= extra
    min_max[1,:] += extra

    colors = []
    # First, draw the borders of all of the boxes
    for box in mesh.boxes:
        colors.append( draw_box(plot, box, min_max) )

    # Draw the centers second in order to make them on 'top' of the borders
    for box,color in zip(mesh.boxes, colors):
        plot.add("%s center"%(box.id), *[[v] for v in box.center],
                 show_in_legend=False, marker_size=5, symbol='square',
                 color=color)



# ====================================================
#      Testing the code for generating box meshes     
# ====================================================

if __name__ == "__main__":
    plot_box_supports = False
    normalize_points = False
    normalize_response = True
    random_points = False
    tight_plot = True # If 'True' will show evalution outside the
    plot_boxes = False #   bounding box for all data points
    plot_true_error = True
    plot_estimated_error = False
    plot_plus_minus_error = False

    dim = 2
    size = 1
    func = linear
    num_points = 8 # Will be to the "dim" power (num_points^dim)
    plot_points = 1000
    plot_range = [[-0.2,1.2]]*dim # Slightly over-sized normal cooridnates

    # if func == quadratic: size += 1/3 # (2/3) / (1/2)
    # elif func == cubic:   size += 2/4 # (3/4) / (1/2)

    #      Various testing functions     
    # ===================================
    # fun = lambda x: 100 * np.sum(x)**2
    # fun = lambda x: 1*x[0]
    # fun = lambda x: (x[0]-num_points/3)*(x[1]-num_points/2)**2
    fun = lambda x: np.cos(x[0]) * (np.sin(x[1]) if len(x) > 1 else 1.0)

    if random_points:
        # Generate testing points randomly spaced
        points = np.random.random(size=(num_points**dim,dim)) * (num_points-1)
    else:
        # Generate testing points in a grid pattern
        points = np.meshgrid(*[range(num_points) for i in range(dim)])
        points = np.array([p.flatten() for p in points]).T

    print(points)

    # Sort the points by y and then x (so they are more intuitively numbered)
    points = points[points[:,1].argsort()]
    points = points[points[:,0].argsort()]
    # Calculate the associated response values
    values = np.array([[fun(pt) for pt in points]]).T

    if normalize_points:
        # Normalize the points themselves
        max_pt_val = float(np.max(points))
        min_pt_val = float(np.min(points))
        points = (points - min_pt_val) / (max_pt_val - min_pt_val)
    else:
        min_val = np.min(points, axis=0)
        max_val = np.max(points, axis=0)
        plot_range = [[ plot_range[i][0] * (max_val[i]-min_val[i]) + min_val[i],
                        plot_range[i][1] * (max_val[i]-min_val[i]) + min_val[i] ]
                      for i in range(dim)]
    if normalize_response:
        # Normalize the response values
        max_resp_val = float(np.max(values))
        min_resp_val = float(np.min(values))
        values = (values - min_resp_val) / (max_resp_val - min_resp_val) + 1.0

    if tight_plot:
        plot_range = [[min(points[:,i]), max(points[:,i])]
                      for i in range(dim)]

    print("Constructing mesh...")
    surf = MaxBoxMesh(func, size)
    surf.fit(points, values)
    print("Creating visualization...")
    surf_p = Plot()
    print(" adding control points...")
    surf_p.add("Control Points", *(np.concatenate((points,values),axis=1).T))
    print(" adding mesh surface...")
    surf_p.add_func("Box Mesh Surface", surf, *plot_range,
                    use_gradient=True, plot_points=plot_points)

    #      Box Mesh Error (Lipschitz) Surfaces     
    # =============================================
    if plot_estimated_error:
        print(" adding estimated relative error surface...")
        surf_p.add_func("Estimated Error", lambda x: sum(surf.error(x)) / abs(surf(x)),
                        *plot_range, plot_points=plot_points)

    if plot_plus_minus_error:
        print(" adding +/- error surface...")
        surf_p.add_func("Min Val", lambda x: surf(x) - surf.error(x),
                        *plot_range, plot_points=plot_points,
                        opacity=0.6)
        surf_p.add_func("Max Val", lambda x: surf(x) + surf.error(x)[1],
                        *plot_range, plot_points=plot_points,
                        opacity=0.6) 

    #      Plotting the True Error     
    # =================================
    if plot_true_error:
        print(" adding true relative error surface...")
        # Create an error function (compensating for normalization)
        def error(x):
            guess = surf(x)
            if normalize_points:
                x = x * (max_pt_val - min_pt_val) + min_pt_val
            resp = fun(x)
            if normalize_response:
               resp = (resp - min_resp_val) / (max_resp_val - min_resp_val) + 1.0
            return abs(guess - resp) / resp

        surf_p.add_func("True Error", lambda x: error(x),
                        *plot_range, plot_points=plot_points)

    #      Plotting the Support Functions     
    # ========================================
    if plot_box_supports:
        print(" adding box supports...")
        for b in surf.boxes:
            surf_p.add_func("%s box"%b.id, b, *plot_range, plot_points=plot_points)

        print(" adding amount of support (*-0.5)")
        #      Support Surface (negated)     
        # ===================================
        surf_p.add_func("Amount of Support", lambda x: -0.5*len([b for b in surf.boxes if b.contains(x)]),
                        *plot_range, use_gradient=True, plot_points=plot_points)

    if plot_boxes:
        print(" creating box plot...")
        boxes_p = Plot()
        draw_boxes(boxes_p, surf)
        print(" making HTML...")
        multiplot([[boxes_p,surf_p]])
    else:
        print(" making HTML...")
        surf_p.plot()
