import random
import plotly
import numpy as np
import colorlover as cl
from scipy.interpolate import PchipInterpolator, splrep, splev
from scipy.spatial import Delaunay, ConvexHull
from scipy.spatial.qhull import QhullError
from scipy.optimize import minimize

PLOT_POINTS = 1000
BRIGHTNESS_RANGE = 0.6

PALATTE = np.array(cl.to_numeric(cl.scales['5']['qual']['Set2']))
# PALATTE = np.array(cl.to_numeric(cl.scales['5']['qual']['Set3']))[::-1]
PALATTE = PALATTE**2
PALATTE = PALATTE / np.max(PALATTE) * 255
# Re-order the palatte so that the colors appear better
PALATTE = np.concatenate((PALATTE[1:], [PALATTE[0]]))
MIN_COLORS = 40

DEFAULT_GRADIENT = np.array(cl.to_numeric(cl.scales['11']['div']['Spectral']))[::-1]

# Expand the palatte
random.seed(0)
palatte_size = len(PALATTE)
for i in range(MIN_COLORS - palatte_size):
    # Create lots of extra colors
    c = np.array([random.choice(PALATTE[:palatte_size,0]),
                  random.choice(PALATTE[:palatte_size,1]),
                  random.choice(PALATTE[:palatte_size,2])])
    # Add this new random color to the palatte
    PALATTE = np.concatenate( (PALATTE, [c]), axis=0 )
# Re-seed the random number generator so that it is not tainted
random.seed()

# PALATTE = np.array([
#     [200, 60, 60],
#     [50, 150,200],
#     [50, 160, 80],
#     [170, 70,250],
#     [220,100, 50],
#     [120,120,120],
# ])

#      Coloring Data     
# =======================

# Given some data, color the data according to a palatte with interpolation
def color_data(values, palatte=DEFAULT_GRADIENT, opacity=1.0):
    no_none = [v for v in values if type(v) != type(None)]
    shift = min(no_none)
    scale = (max(no_none) - shift) * 1.11
    def color(value):
        if value == None: return None
        index = len(palatte) * (value-shift) / scale
        lower = int(index)
        upper = lower + 1
        if (lower > len(palatte)-1): lower = len(palatte)-1
        if (upper > len(palatte)-1): upper = len(palatte)-1
        index -= lower
        c = tuple(palatte[lower]*(1-index) + palatte[upper]*(index))
        return 'rgba(%i,%i,%i,%f)'%(c+(opacity,))
    return list(map(color, values))

#      Functional Representations of Data     
# ============================================


# Given three points, this will solve the equation for the quadratic
# function which intercepts all 3
def solve_quadratic(x, y):
    if len(x) != len(y): raise(Exception("X and Y must be the same length."))
    if len(x) != 3: raise(Exception("Exactly 3 (x,y) coordinates must be given."))
    x1, x2, x3 = x
    y1, y2, y3 = y
    a = -((-x2 * y1 + x3 * y1 + x1 * y2 - x3 * y2 - x1 * y3 + x2 * y3)/((-x1 + x2) * (x2 - x3) * (-x1 + x3)))
    b = -(( x2**2 * y1 - x3**2 * y1 - x1**2 * y2 + x3**2 * y2 + x1**2 * y3 - x2**2 * y3)/((x1 - x2) * (x1 - x3) * (x2 - x3)))
    c = -((-x2**2 * x3 * y1 + x2 * x3**2 * y1 + x1**2 * x3 * y2 - x1 * x3**2 * y2 - x1**2 * x2 * y3 + x1 * x2**2 * y3)/((x1 - x2) * (x1 - x3) * (x2 - x3)))
    return (a,b,c)

# Returns the set of mode points given the CDF second difference list,
# mode points are sorted in order of magnitude of relative differences
def mode_points(second_diff, thresh_perc=0.0):
    # Collect together the maximum absolute values in each
    # group of signs (find local peaks)
    sign_maxes = []
    curr_max = second_diff[0]
    thresh = np.percentile(second_diff[:,1], thresh_perc)
    for val in second_diff:
        if (abs(val[1]) > thresh) and (val[1] * curr_max[1] < 0):
            sign_maxes.append(curr_max)
            curr_max = val
        elif abs(val[1]) > abs(curr_max[1]):
            curr_max = val
    sign_maxes = np.array(sign_maxes)
    # Record the magnitude of the difference between the peaks
    # as well as the location in the middle of the two peaks    
    sign_changes = []
    for i in range(1,len(sign_maxes)):
        x_val = (sign_maxes[i][0] + sign_maxes[i-1][0]) / 2
        y_val = abs(sign_maxes[i][1] - sign_maxes[i-1][1])
        sign_changes.append([x_val, y_val])
    # Sort the valley
    sign_changes.sort(key=lambda i: -i[1])
    sign_changes = np.array(sign_changes)
    shift = min(sign_changes[:,1])
    scale = max(sign_changes[:,1]) - shift
    sign_changes[:,1] = (sign_changes[:,1] - shift) / scale
    return sign_changes

# CDF Second Difference
def cdf_second_diff(data):
    # Sort the data and get the min and max
    data = list(data); data.sort()
    data = np.array(data)
    min_pt = data[0]
    max_pt = data[-1]
    # Initialize a holder for all CDF points
    data_pts = []
    # Add all of the CDF points for the data set
    for i,val in enumerate(data):
        if ((i > 0) and (val == data[i-1])): continue
        data_pts.append( [val, i/(len(data) - 1)] )
    # Add the 100 percentile point if it was not added already
    if data_pts[-1][1] != 1.0: data_pts[-1][1] = 1.0
    # Convert it all to numpy format for ease of processing
    data_pts = np.array(data_pts)
    second_diff_pts = []
    for i in range(1,len(data_pts)-1):
        a,_,_ = solve_quadratic(data_pts[i-1:i+2,0],
                                data_pts[i-1:i+2,1])
        second_deriv = 2*a
        second_diff_pts.append( [data_pts[i,0], second_deriv] )
    # # Sort the second_diff points by the magnitude of the second derivative
    # second_diff_pts.sort(key=lambda pt: -abs(pt[1]))
    return np.array(second_diff_pts)
        
        
# Linearly interpolate to guess y values for x between provided data
def linear_fit_func(x_points, y_points, min_max=None):
    # Generate the 'range' of the function if it is not provided
    if type(min_max) == type(None):
        min_max = (min(x_points), max(x_points))

    # Fit the data with degree 1 splines (linear interpolation)
    fit = splrep(x_points, y_points, k=1)

    # Generate the function to return (gives the 'range' of
    # application when called without x values)
    def func(x_val=None, fit=fit, min_max=min_max):
        if type(x_val) == type(None):
            return min_max
        else:
            return splev(x_val, fit, ext=3)
    def deriv(x_val, fit=fit):
        return splev(x_val, fit, der=1, ext=3)

    # Make the attribute "derivative" return a function for the derivative
    setattr(func, "derivative", lambda: deriv)

    # Return the linear interpolation function
    return func

# Subsample a set of data "subsamples" times with each sample of size
# "subsample_size" and then generate a list of fit-functions for each
# of the percentiles in "percentiles". Return list of fit functions.
def percentile_funcs(x_points, y_points, percentiles, min_max):
    # Generate the points for the percentile CDF's
    perc_points = []    
    for p in percentiles:
        perc_points.append(
            [np.percentile(y_points[:,i], p) for i in range(len(x_points))]
        )
    perc_points = np.array(perc_points)

    # Generate the CDF functions for each percentile
    funcs = []
    for pts in perc_points:
        # Generate a linear fit function over each set of percentile points
        func = linear_fit_func(x_points, pts)
        funcs.append(func)

    # Return the fit functions for each of the percentiles
    return funcs

# Returns the CDF function value for any given x-value
def cdf_fit_func(data, linear=False):
    # Sort the data and get the min and max
    data.sort()
    min_pt = data[0]
    max_pt = data[-1]

    # Initialize a holder for all CDF points
    data_pts = []

    # Add all of the CDF points for the data set
    for i,val in enumerate(data):
        if ((i > 0) and (val == data[i-1])): continue
        data_pts.append( [val, i/(len(data) - 1)] )
    # Add the 100 percentile point if it was not added already
    if data_pts[-1][1] != 1.0: data_pts[-1][1] = 1.0

    # Convert it all to numpy format for ease of processing
    data_pts = np.array(data_pts)

    # Generate a fit for the data points
    if linear:
        fit = splrep(data_pts[:,0], data_pts[:,1], k=1)
        fit = lambda x_val, fit=fit: splev(x_val, fit)
        fit.derivative = lambda d: (lambda x_val: splev(x_val, fit, der=d))
    else:
        fit = PchipInterpolator(data_pts[:,0], data_pts[:,1])
    
    # Generate a function that computes this CDF's points
    def cdf_func(x_val=None):
        if type(x_val) == type(None):
            return (min_pt, max_pt)        
        else:
            y_val = fit(x_val)
            if (type(x_val) == type(np.array([]))):
                y_val = np.where(x_val < min_pt, 0.0, y_val)
                y_val = np.where(x_val > max_pt, 1.0, y_val)
            return y_val
    def deriv(x_val):
        val = fit.derivative(1)(x_val)
        val = np.where(x_val < min_pt, 0, val)
        val = np.where(x_val > max_pt, 0, val)
        return val

    setattr(cdf_func, "derivative", lambda: deriv)

    # Return the custom function for this set of points
    return cdf_func

# Returns the PDF function for the data
def pdf_fit_func(data=None):
    # Take the first derivative of the CDF function to get the PDF
    return cdf_fit_func(data).derivative(1)

# Calculate the maximum difference between two CDF functions (two sample)
def ks_diff(test_func, true_func):
    # Cycle through the functions to find the min and max of all ranges
    min_pt, max_pt = true_func()

    diff_func = lambda x: -abs(test_func(x) - true_func(x))
    # Use scipy minimize to find the greatest difference between the functions
    # sol = minimize(diff_func, [(max_pt - min_pt) / 2],
    #                method='L-BFGS-B', bounds=[(min_pt,max_pt)]).x

    # Generate a large set of x-points (for a smooth plot)
    x_points = np.linspace(min_pt, max_pt, 1000)
    diff = abs(test_func(x_points) - true_func(x_points))

    greatest_diff = -float('inf')
    while (diff[np.argmax(diff)] > greatest_diff):
        lower = np.argmax(diff) - 1
        upper = np.argmax(diff) + 1
        min_pt = x_points[max(lower, 0)] - (
            1 if lower < 0 else 0)
        max_pt = x_points[min(upper,len(x_points)-1)] + (
            1 if upper >= len(x_points) else 0)
        x_points = np.linspace(min_pt, max_pt, 1000)
        diff = abs(test_func(x_points) - true_func(x_points))
        greatest_diff = max(max(diff), greatest_diff)

    return greatest_diff


#      Class for building a plotly plot     
# ==========================================

class Plot:
    def __init__(self, title="", x_title="x", y_title="y",
                 z_title="z", x_range=None, y_range=None,
                 z_range=None, mode="markers", palatte=PALATTE):
        self.title = title
        self.x_title = x_title
        self.y_title = y_title
        self.z_title = z_title
        self.x_min_max = [float('inf'), -float('inf')]
        self.y_min_max = [float('inf'), -float('inf')]
        self.z_min_max = [float('inf'), -float('inf')]
        self.x_range = x_range
        self.y_range = y_range
        self.z_range = z_range
        # Specific booleans for tracking internal state
        self.is_3d = False
        self.to_reverse = []
        # Data for tracking default plot settings
        self.color_num = -1
        self.data = []
        self.mode = mode
        self.palatte = palatte
        self.palatte_size = len(palatte)

    # Get a color from a palatte
    def color(self, number=None, brightness=1.0, alpha=1.0, color=None):
        if color == None:
            if (number == None): raise(Exception("ERROR: Neither color nor number have been provided."))
            if (number < len(self.palatte)):
                # If we have fewer entries than the palette size
                c = self.palatte[number]
            else:
                # Otherwise we have to create a new palette entry
                c = np.array([random.choice(self.palatte[:self.palatte_size,0]),
                              random.choice(self.palatte[:self.palatte_size,1]),
                              random.choice(self.palatte[:self.palatte_size,2])])
                # Add this new random color to the palatte
                self.palatte = np.concatenate( (self.palatte, [c]), axis=0 )
        else:
            # Get the color as a list of numbers
            c = color[color.index('(')+1:color.index(')')].split(',')
            # Make sure the color only has [red, green, blue, alpha]
            c = np.array(list(map(float,c)))[:3]

        # Apply the brightness to the color
        c = c*brightness
        c = np.where(c > 255, 255, c)
        c = np.where(c < 0, 0, c)

        # Return the color as a plotly color string
        return 'rgba(%i,%i,%i,%f)'%(tuple(c)+(alpha,))

    # Return the face color of a simplex given the simplex, data z
    # values, and data index
    def simp_color(self, simp, z, color_ind, opacity, colors=None):
        no_none = [v for v in z if type(v) != type(None)]
        shift = min(no_none)
        width = max(no_none) - shift
        has_none = type(None) in (type(v) for v in z[simp])
        if (width != 0) and (not has_none): 
            # If colors were provided, then average them to produce out color
            if type(colors) != type(None):
                # Get the color of each node in the simplex as a numpy array
                colors = [colors[i] for i in simp]
                colors = [c[c.index('(')+1:c.index(')')].split(',')
                          for c in colors]
                colors = np.array([list(map(float,c)) for c in colors])
                if colors.shape[1] != 4:
                    colors = np.concatenate((
                        colors,np.ones(shape=(colors.shape[0],1))),
                                            axis=1)
                # return the average color of points in the simplex
                return 'rgba(%f,%f,%f,%f)'%tuple(np.sum(colors,axis=0) / len(simp))
            else:
                simp_avg = sum(z[simp]) / len(simp)
                brightness = (1.0-BRIGHTNESS_RANGE/2) + ((simp_avg - shift) / width) * BRIGHTNESS_RANGE
        else:
            brightness = 1.0
        return self.color(color_ind, brightness, opacity)


    # "func" returns True or False if a point is in the region
    # ONLY WORKS FOR 2D DATA
    def add_region(self, name, func, min_max_x=None, min_max_y=None,
                   plot_points=PLOT_POINTS, **kwargs):
        if self.is_3d: raise(Exception("ERROR: Regions only work for 2D plots."))
        if type(min_max_x) == type(None):
            min_max_x = self.x_min_max.copy()
        if type(min_max_y) == type(None):
            min_max_y = self.y_min_max.copy()
        if max(map(abs,min_max_x+min_max_y)) == float('inf'):
            raise(Exception("ERROR: Invalid x an y range."))
        # Round up the number of plot points per axis
        plot_points = int(plot_points**(0.5) + 0.5)
        # Calculate the mesh grid of x and y values
        x_vals = (np.linspace(*min_max_x, num=plot_points),)
        y_vals = (np.linspace(*min_max_y, num=plot_points),)
        x,y = np.meshgrid(x_vals, y_vals)
        test_pts = np.vstack((x.flatten(), y.flatten())).T
        region_pts = np.array([pt for pt in test_pts if func(pt)])
        try:
            hull_pts = region_pts[ConvexHull(region_pts).vertices]
            self.add(name, hull_pts[:,0], hull_pts[:,1], mode='lines',
                     fill='toself', line_width=0, opacity=0.1, **kwargs)
        except:
            print("ERROR: Could not add region.")


    # Adds a function to the plot given minimum and maximum values
    def add_func(self, name, func, min_max_x, min_max_y=[],
                 mode=None, plot_type=None, plot_points=PLOT_POINTS,
                 grid_lines=True, use_gradient=False, **kwargs):
        if len(min_max_y) > 0: self.is_3d = True
        # If we have two control axes, square root the plot points
        if self.is_3d:
            plot_points = int(plot_points**(0.5) + 0.5)
            # If no y was provided, set it to default value
            if len(min_max_y) == 0: min_max_y = [0.0,0.0]
            if mode == None: plot_type = 'surface'
        else:
            if mode == None: mode = 'lines'

        # Generate the input points
        x_vals = (np.linspace(*min_max_x, num=plot_points),)
        if self.is_3d:
            x_vals += (np.linspace(*min_max_y, num=plot_points),)
        x_vals = tuple(x.flatten() for x in np.meshgrid(*x_vals))

        # Get the response values
        try:
            response = np.array([func(x) for x in np.vstack(x_vals).T]).flatten()
        except SystemExit: exit()
        except:
            # Provide a useful error message if the response values
            # could not be computed correctly
            try:
                sample = func((np.vstack(x_vals).T)[0])
                raise(Exception("Error in return value of provided function. Expected number, got '%s'"%(type(sample))))
            except:
                raise(Exception("Error computing the provided function."))

        # Make a nice pretty gradient of color
        if ('color' not in kwargs) and use_gradient:
            colors = color_data(response)
            kwargs['marker_colors'] = colors

        # Call the standard plot function
        self.add(name, *x_vals, response, mode=mode,
                 plot_type=plot_type, **kwargs)
            
        # If this is a 3D surface plot, add grid lines
        if (self.is_3d and plot_type == 'surface') and grid_lines:
            opacity = kwargs.get("opacity",1.0)
            line_color = 'rgb(0,0,0)'
            y_vals = np.array(list(range(plot_points))) * plot_points
            for row in range(plot_points):
                x = x_vals[0][row*plot_points:(row+1)*plot_points]
                y = x_vals[1][row*plot_points:(row+1)*plot_points]
                z = response[row*plot_points:(row+1)*plot_points]
                self.add("", x,y,z, show_in_legend=False,
                         group=name+"lines", mode='lines',
                         color=line_color, line_width=opacity, opacity=opacity)

                indices = [row + i*plot_points for i in range(plot_points)]
                x = x_vals[0][indices]
                y = x_vals[1][indices]
                z = response[indices]
                self.add("", x,y,z, show_in_legend=False,
                         group=name+" (lines)", mode='lines',
                         line_width=opacity, opacity=opacity,
                         color=line_color)


    # Function for adding a new series to the plot
    def add(self, name, x_values=None, y_values=None, z_values=None,
            text=None, plot_type=None, group=None,
            show_in_legend=True, fill=None, line_width=None,
            line_color=None, fill_color=None, mode=None, dash=None,
            symbol='circle', color=None, opacity=1.0,
            fill_opacity=0.6, marker_size=None, marker_colors=None,
            shade=True, marker_line_width=0,
            marker_line_color='rgba(50,50,50,0.8)', **kwargs):

        # Convert the x, y (and z) values into numpy arrays and
        # store 'values' for creating marker colors based on magnitude
        if type(x_values) != type(None):
            x_values = np.array(x_values)
            values = x_values
            no_none = [v for v in x_values if type(v) != type(None)]
            if len(no_none) != 0:
                self.x_min_max = [min(min(no_none), self.x_min_max[0]),
                                  max(max(no_none), self.x_min_max[1])]
        if type(y_values) != type(None):
            y_values = np.array(y_values)
            values = y_values
            no_none = [v for v in y_values if type(v) != type(None)]
            if len(no_none) != 0:
                self.y_min_max = [min(min(no_none), self.y_min_max[0]),
                                  max(max(no_none), self.y_min_max[1])]
        if type(z_values) != type(None):
            self.is_3d = True
            z_values = np.array(z_values)
            values = z_values
            no_none = [v for v in z_values if type(v) != type(None)]
            if len(no_none) != 0:
                self.z_min_max = [min(min(no_none), self.z_min_max[0]),
                                  max(max(no_none), self.z_min_max[1])]

        # Define z-values if none were given and we need them, and plot type
        if self.is_3d:
            if plot_type == None:
                plot_type = 'scatter3d'
            if type(z_values) == type(None):
                z_values = np.zeros(len(x_values))
            if text == None:
                text = ["%s: %s<br>%s: %s<br>%s: %s"%(
                    self.x_title,x, self.y_title,y, self.z_title,z)
                        for (x,y,z) in zip(x_values,y_values,z_values)]
        else:
            if plot_type == None:
                plot_type = 'scatter'

        # Process mode
        if type(mode) == type(None):
            mode = self.mode
        # Set the color if none was provided
        if color == None:
            self.color_num += 1
            color = self.color(self.color_num, alpha=opacity)
        else:
            shade = False
        if line_color == None:
            line_color = color
        if fill_color == None:
            fill_color = color
        if not marker_colors:
            if shade:
                marker_colors = []
                no_none = [v for v in values if v != None]
                shift = min(no_none)
                scale = max(no_none) - shift
                if scale == 0: scale = 1.0
                for v in values:
                    brightness = (1.0-BRIGHTNESS_RANGE/2) + ((v - shift) / scale) * BRIGHTNESS_RANGE
                    marker_colors.append( self.color(color=color, brightness=brightness, alpha=opacity) )
            else:
                marker_colors = [color]*len(values)

        # Special plotly failure mode, need to reverse data for
        # 'tonext' to actually mean 'next' instead of 'previous'
        self.to_reverse.append(
            (type(fill) == type("")) and ("tonext" in fill))

        # Now add the data to local storage
        self.data.append(dict(
            type = plot_type,
            name = name,
            x = x_values,
            y = y_values,
            z = z_values,
            hoverinfo = 'text',
            text = text,
            marker = dict(
                # Generate colors based on point magnitude
                color = color if ("lines" in mode) else marker_colors,
                # color = marker_colors,
                size = marker_size,
                opacity = opacity,
                symbol = symbol,
                line = dict(
                    width = marker_line_width,
                    color = marker_line_color
                )),
            line = dict(
                width = line_width,
                color = line_color,
                dash = dash
            ),
            mode = mode,
            fill = fill,
            fillcolor = fill_color,
            legendgroup = group,    
            showlegend = show_in_legend
        ))

        # Update the newly created dictionary with any custom user settings
        self.data[-1].update(kwargs)

    # Prepares all the data sets to be plotted in whatever dimension
    # is highest (2 or 3). Creates 3D meshes for all surfaces.
    def clean_data(self):
        # Remove all references to 3D layout if this is a 2D plot
        if not self.is_3d:
            # 2D PLOT SETUP
            for d in self.data:
                d.pop('z','')
                d.pop('hoverinfo','')
                d.pop('text','')
                # Special case for plotting histograms
                if d['type'] == 'histogram':
                    if type(d.get('y','')) == type(None):
                        d.pop('y','')
                    if type(d.get('x','')) == type(None):
                        d.pop('x','')
                    d.pop('line','')
                    d.pop('mode','')
                    d.pop('fill','')
                    d['opacity'] = d['marker'].pop('opacity','')
                    d['marker'].pop('symbol','')
                    d['marker'].pop('size','')
                    d['marker']['color'] = d.pop('fillcolor','')
        else:
            # 3D PLOT SETUP
            for ind,d in enumerate(self.data):
                # Add z values to all scatters that may have been added
                if d['type'] == 'scatter':
                    d['z'] = np.zeros(len(d['x']))
                    d['type'] = 'scatter3d'                            
                    if d['marker']['size'] == None:
                        d['marker']['size'] = 5
                # Convert fill and / or lines into surfaces
                conv_2d = (not self.is_3d) and ('lines' in d['mode'])
                if (d.get('fill','') == 'toself') or conv_2d:
                    print("CONVERTING 2D")
                    d['type'] = 'surface'
                    # Get the opacity of the surface
                    if d.get('fill','') != None:
                        d['opacity'] = float(d['fillcolor'].split(',')[-1].strip(')'))
                    else:
                        d['opacity'] = float(d['line']['color'].split(',')[-1].strip(')'))
                # If user wants a surface, construct one! 
                # (plotly default surfaces are not very cooperative)
                if ('surface' in d['type']):
                    points_2D = np.vstack([d['x'], d['y']]).T
                    try:
                        simps = Delaunay(points_2D).simplices
                        d['type'] = 'mesh3d'
                        d['i'] = simps[:,0]
                        d['j'] = simps[:,1]
                        d['k'] = simps[:,2]
                        # Generate face colors with average simplex z-coordinate
                        d['facecolor'] = list(map(
                            lambda simp: self.simp_color(
                                simp, d['z'], ind,
                                d['marker']['opacity'],
                                d['marker']['color']), simps
                        ))
                        if 'opacity' not in d:
                            d['opacity'] = d['marker']['opacity']
                        # Currently there's a bug in plotly that doesn't allow text to be set manuall, this is the only way to get hover info
                        d['hoverinfo'] = 'all'
                        d.pop('marker','')
                        d.pop('mode','')
                        d.pop('text','')
                    except QhullError:
                        d['type'] = 'scatter3d'
                        if 'mode' not in d:
                            d['mode'] = 'lines'
                # Pop out the unnecessary attributes for 3D plots
                d.pop('fill','')
                d.pop('fillcolor','')
                if 'line' not in d.get('mode',''):
                    d.pop('line','')

    # Manage plotly reverse order bug (only happens with "tonext_")
    def reorder_data(self):
        start = end = None
        # Cycle through the elements of data
        for i,tr in enumerate(self.to_reverse):
            if (tr and start==None):
                start = i
            if (not tr and start!=None):
                end = i+1
                # Reverse that group of plot series
                self.data[start:end] = self.data[start:end][::-1]
                start = end = None
        # Reverse the final group when self.to_reverse[-1] == True
        if (start!=None):
            end = len(self.data)
            self.data[start:end] = self.data[start:end][::-1]

        self.to_reverse = [False] * len(self.data)

    # Alias to the "plot" function
    def draw(self, *args, **kwargs):
        return self.plot(*args, **kwargs)

    # Function for plotting the data that has been added simply,
    # handles most things automatically, particularly whether or not
    # the user is plotting in 3D
    def plot(self, title=None, x_range=None, y_range=None,
             z_range=None, show=True, file_name="temp-plot.html",
             fixed=True, show_legend=True, html=True, layout={},
             scene_settings={}, **kwargs):
        # Update title, and all plot axis ranges
        if title == None:
            title = self.title
        if (fixed and x_range == None and
            max(map(abs,self.x_min_max)) != float('inf')):
            x_width = self.x_min_max[1] - self.x_min_max[0]
            x_range = [self.x_min_max[0] - 0.05*x_width,
                       self.x_min_max[1] + 0.05*x_width]
        if (fixed and y_range == None and
            max(map(abs,self.y_min_max)) != float('inf')):
            y_width = self.y_min_max[1] - self.y_min_max[0]
            y_range = [self.y_min_max[0] - 0.05*y_width,
                       self.y_min_max[1] + 0.05*y_width]
        if (fixed and z_range == None and
            max(map(abs,self.z_min_max)) != float('inf')):
            z_width = self.z_min_max[1] - self.z_min_max[0]
            z_range = [self.z_min_max[0] - 0.05*z_width,
                       self.z_min_max[1] + 0.05*z_width]
        # Generate the layout (titles and legend)
        plot_layout = dict(
            title = title,
            showlegend = show_legend,
        )
        # Setup for the axes of the plot
        scene = dict(
            xaxis = dict(title = self.x_title, range=x_range),
            yaxis = dict(title = self.y_title, range=y_range),
            zaxis = dict(title = self.z_title, range=z_range),
        )
        # Setup the plot layout (different for 2D and 3D plots)
        if not self.is_3d:
            plot_layout.update(scene)
            plot_layout.pop('zaxis')
        else:
            scene['aspectmode'] = 'cube'
            scene['camera'] = dict(
                up=dict(x=0, y=0, z=1),
                center=dict(x=0, y=0, z=0),
                eye=dict(x=-1.0, y=-2.0, z=0.7)
            )
            scene.update(scene_settings)
            plot_layout['scene'] = scene

        # Make sure all the data entries are prepared to be plotted
        self.clean_data()
        # Update the plot layout with any specific user settings given
        plot_layout.update(layout)         
        data = self.data

        # Manage plotly reverse order bug (only happens with "tonext_")
        self.reorder_data()
        # Generate the figure
        fig = plotly.graph_objs.Figure(data=data, layout=plot_layout)

        # Create the html file and show in browser if appropriate
        if html: create_html(fig, file_name, show, **kwargs)
        # Return the figure
        return fig

        
#      Functions for manipulation produces plots     
# ===================================================

# Generates the HTML file and fixes some javascript so that the
# plot does not have unwanted buttons and links
def create_html(fig=None, file_name="temp-plot.html", show=True, **kwargs):
    # Generate the plot offline 
    plotly.offline.plot(fig, auto_open=show, filename=file_name,
                        show_link=False, **kwargs)
    # Remove unnecessary modebar buttons and the plotly logo link
    with open(file_name) as f:
        file_string = f.read()
    file_string = file_string.replace(
        'displaylogo:!0', 'displaylogo:!1')
    file_string = file_string.replace(
        'modeBarButtonsToRemove:[]',
        'modeBarButtonsToRemove:["sendDataToCloud", "select2d", "lasso2d"]')
    file_string += "\n\n"
    with open(file_name, "w") as f:
        f.write(file_string)
    return file_name


# Make multiple plots fit onto one browser window, options for sharing
# axes as well for naming purposes. Mixed plot types allowed too!
# Supports different number of columns per row, but no other mismatch
def multiplot(plots, x_domains=None, y_domains=None, show=True, specs=None,
              shared_y=False, shared_x=False, gap=0.05, **kwargs):
    # Make sure the plots array is 2D
    try:    plots[0][0]
    except: plots = [plots]
    # Convert given plots into figures (if figures were not given
    for r in plots:
        for c in range(len(r)):
            if type(r[c]) == Plot:
                r[c] = r[c].plot(html=False, show=False)

    # Count the number of rows and columns
    rows = len(plots)
    cols = [len(r) for r in plots]
    max_cols = max(c for c in cols)
    # Generate/Process the specs
    if type(specs) != type(None):
        try:    specs[0][0]
        except: specs = [specs]
    else:
        specs = [[None]*max_cols for r in range(rows)]
        for r,row in enumerate(plots):
            for c,plot in enumerate(row):
                if type(plot) == type(None): continue
                sample_data = plots[r][c]['data'][0]
                specs[r][c] = {"is_3d": 'z' in sample_data}
    # Generate the x and y domains if they are not provided by the user
    if x_domains == None:                
        x_domains = []
        for r in range(rows):
            width = (1 - (cols[r]-1)*gap) / cols[r]
            x_domains.append(
                [[c*(width+gap), c*(width+gap) + width]
                 for c in range(cols[r])])
    if y_domains == None:                
        height = (1 - (rows-1)*gap) / rows
        y_domains = [[r*(height+gap), r*(height+gap) + height]
                     for r in range(rows)]
    # Identify the number of dimensions provided in x an y domains, if
    # too few, then make sure it is the same shape as the plots
    try: x_domains[0][0][0]
    except TypeError:
        x_domains = [x_domains for r in range(rows)]
    try: y_domains[0][0][0]
    except TypeError:
        y_domains = [[y_domains[r]]*cols[r] for r in range(rows)]
            
    # Fix y-domains so that they are specified from bottom to top
    flipped_y = []
    gap = y_domains[1][0][0] - y_domains[0][0][1] if len(y_domains) > 1 else 0
    for r in range(rows):
        start = 0.0 if r == 0 else flipped_y[-1][1] + gap
        width = y_domains[rows-r-1][0][1] - y_domains[rows-r-1][0][0]
        flipped_y.append([start, start+width])
    y_domains = [[flipped_y[r]]*cols[len(cols)-1-r] for r in range(rows)][::-1]

    # Generate the holder for the multiplot
    fig = plotly.tools.make_subplots(rows=rows, cols=max_cols,
                                     specs=specs,
                                     shared_yaxes=shared_y,
                                     shared_xaxes=shared_x)
    # Generate the multi plot!
    counter_2d = 0
    counter_3d = 0
    for r,row in enumerate(plots):
        for c,plot in enumerate(row):
            # Allows for empty spaces
            if type(plot) == type(None): continue
            # Otherwise, continue assuming we have a figure!
            for d in plot['data']:
                fig.append_trace(d, r+1, c+1)

            if specs[r][c]['is_3d']:
                counter_3d += 1
                scene_name = 'scene' + str(counter_3d)
                fig['layout'][scene_name].update(plot['layout']['scene'])
                fig['layout'][scene_name]['domain']['x'] = x_domains[r][c]
                fig['layout'][scene_name]['domain']['y'] = y_domains[r][c]

            else:
                counter_2d += 1
                x_name = 'xaxis'+str(counter_2d)
                y_name = 'yaxis'+str(counter_2d)
                # For shared axes, only add the first entry of column or row
                # Update the domains as specified by the user
                if (not shared_x) or (r == 0):
                    fig['layout'][x_name].update(plot['layout'].pop('xaxis'))
                    fig['layout'][x_name]['domain'] = x_domains[r][c]
                if (not shared_y) or (c == 0):
                    fig['layout'][y_name].update(plot['layout'].pop('yaxis'))
                    fig['layout'][y_name]['domain'] = y_domains[r][c]
                # Ensure that no axis layouts make it into the plot that shouldn't
                plot['layout'].pop('xaxis','')
                plot['layout'].pop('yaxis','')
            fig['layout'].update(plot['layout'])
            # Remove the 'scene' if there is one left over
            if specs[r][c]['is_3d']: fig['layout'].pop('scene','')

    # Create the html plot if the user wants that (pass extra arguments)
    if show: create_html(fig, **kwargs)
    # Return the figure to be plotted
    return fig



# Given a plot, this will add the 'percentiles' cloud
def plot_percentiles(plot, name, funcs, percentiles, min_max_x,
                     center_color = np.array([255,70,0,0.6]),
                     outer_color = np.array([255,150,50,0.3])):
    group_id = str(np.random.randint(1000000))
    percentiles.sort() # Ensure they are sorted
    gaps = [percentiles[i] - percentiles[i-1] for i in range(1,len(percentiles))]
    textfont = dict(color='rgb(40,40,40)')
    for f,p,g in zip(funcs[:-1], percentiles[:-1], gaps):
        ratio = abs((50 - p - g/2)/50)
        color = center_color * abs(1.0 - ratio) + \
                outer_color  * ratio
        color = 'rgba(%i,%i,%i,%f)'%tuple(color)
        plot.add_func("%ith percentile"%p, f, min_max_x, fill='tonexty',
                      fill_color=color, mode='none', group=group_id,
                      show_in_legend=False, textfont=textfont)

    # Add the last function (bounding the set)
    plot.add_func("%ith percentile"%percentiles[-1], funcs[-1],
                  min_max_x, mode='lines', line_width=0,
                  group=group_id, show_in_legend=False, color=color,
                  textfont=textfont)

    # Add a master legend entry
    x_val = (min_max_x[1] - min_max_x[0])/2
    y_val = funcs[len(funcs)//2](x_val)
    plot.add(name + " Percentiles", [x_val], [y_val], mode='lines',
             color='rgba(%i,%i,%i,%f)'%tuple(center_color),
             group=group_id)

    return plot



if __name__ == "__main__":

    # # Testing out the mathematical routines
    # num_modes = 5
    # top_mode_color = "rgba(20,20,20,0.6)"
    # rest_mode_color = "rgba(210,210,210,0.4)"
    # size = 200
    # values = np.vstack((np.random.normal(0.4, 0.05, size=(size,)),
    #                     np.random.normal(0.0, 0.02, size=(size,)))).flatten()
    # values.sort()

    # raw = Plot("", "Value", "Normalized Scale")
    # cdf = Plot("", "Value", "Percent Less")
    # pdf = Plot("CDF and PDF of random normal 2-mode data", "Value", "Probability of Occurrence")

    # min_max = (min(values), max(values))
    # bin_size = (min_max[1]-min_max[0])/100
    # pdf.add("PDF Histogram", x_values=values, plot_type="histogram",group="1",
    #         opacity=0.7, autobinx=False, histnorm='probability',
    #         xbins=dict(start=min_max[0],end=min_max[1],size=bin_size))
    # cdf.add_func("CDF", cdf_fit_func(values, linear=True), min_max,
    #              mode="markers+lines", marker_size=4, group="1")

    # sd_pts = cdf_second_diff(values)

    # # Add vertical lines for each of the mode locations
    # mode_pts = mode_points(sd_pts, thresh_perc=50.0)
    # # Add the rest of the modes as another series
    # mode_x_vals = []
    # mode_y_vals = []
    # for pt in mode_pts[num_modes:]:
    #     mode_x_vals += [pt[0], pt[0], None]
    #     mode_y_vals += [0,     1,     None]
    # cdf.add("Remaining %i Modes"%(len(mode_pts) - num_modes),
    #         mode_x_vals, mode_y_vals, mode='lines', line_width=1,
    #         color=rest_mode_color, group='modes')
    # # Add the 'num_modes' top modes as a series
    # mode_x_vals = []
    # mode_y_vals = []
    # for pt in mode_pts[:num_modes]:
    #     mode_x_vals += [pt[0], pt[0], None]
    #     mode_y_vals += [0,     1,     None]
    # cdf.add("Top %i Modes"%num_modes,mode_x_vals, mode_y_vals,
    #         mode='lines', line_width=1, color=top_mode_color, group='modes')


    # # Generate the discrete CDF second derivative plot
    # mode_pts = [list(pt) for pt in mode_pts]
    # mode_pts.sort(key=lambda pt: pt[0])
    # mode_pts = np.array(mode_pts)
    # raw.add("Absolute Discrete CDF''", mode_pts[:,0], mode_pts[:,1],
    #         mode='lines', fill='tozeroy', color=top_mode_color,
    #         fill_color=rest_mode_color, line_width=1, group='modes')


    # # Plot all three of these stacked together
    # multiplot([[raw.plot(show=False, html=False)],
    #            [cdf.plot(show=False, html=False)],
    #            [pdf.plot(show=False, html=False)]
    # ], shared_x=True)
    
    # exit()


    # Testing code for the plotting interface
    fun = lambda x: np.sum(x**2) / 10
    # fun = lambda x: x[-1]*x[-2]
    x = np.linspace(-10,10,100)
    y = x**2 / 10    

    # Decide whether to make a 2D, 3D, or both
    selection = np.random.randint(3)
    selection = 2
    
    # How to add a histogram to the plot

    # pdf_plot.add("<name>", y_values=<list of values>, plot_type="histogram",
    #              opacity=0.7, autobiny=False,
    #              histnorm='probability',
    #              ybins=dict(start=<min val>,end=<max val>,size=<bin size>))

    # Plot the 2D
    if selection in [0,2]:
        plot = Plot()        
        plot.add_func("Test Func 2D", fun,[-10,10], opacity=0.5)
        # plot.add_func("Test Func 2D", lambda x: fun(x) - 3,[-10],[10], opacity=0.5)
        # plot.add_func("Test Func 2D", lambda x: fun(x) - 2,[-5],[5], opacity=0.5)
        # plot.add("Test 1", x, y)
        plot.add("V Line", [0,0], [min(y), max(y)], mode="lines+markers")
        plot.add("Square", [-2,-2,2,2], [5,10,10,5],
                 mode="none", fill="toself")
        plot.add("H Line", [-5,5], [1,1], mode="lines+markers", symbol='square')
        plot.add("H Line 2", [-5,5], [2,2], mode="lines")
        plot1 = plot

    # Plot the 3D
    if selection in [1,2]:
        plot = Plot("Title","X Axis", "Y Axis", "Z Axis")
        rand_x = list(range(-5,6,2))
        rand_y = np.random.randint(-3,4,size=6)
        rand_z = np.random.randint(3,8,size=6)
        plot.add("3D Points", rand_x, rand_y, rand_z, mode='lines', fill='tozeroy')
        dens = 5
        x, y = np.meshgrid(np.linspace(-5,5,dens), np.linspace(-5,5,dens))
        x = x.flatten()
        y = y.flatten()
        fun = lambda x: -.3*x[1] + 1/2*x[0] + 1
        z = np.array(list(map(fun, zip(x,y))))

        plot.add("3D Above", x, y, z+1.5, marker_line_width=1)
        plot.add("3D Below", x, y, z-1.5, marker_line_width=1)
        plot.add("3D B Next", x, y, z-1, opacity=0.7, marker_line_width=1)
        plot.add("3D A Next", x, y, z+1, opacity=0.7, marker_line_width=1)
        plot.add_func("3D Surface", fun, [min(x),max(x)], [min(y),max(y)], opacity=0.7)
        plot2 = plot


    fig1 = plot1.plot(show=False)
    fig2 = plot2.plot(show=False)

    specs = [{"is_3d":plot1.is_3d}, {"is_3d":plot2.is_3d}]

    multiplot([fig1, fig2], specs=specs)

