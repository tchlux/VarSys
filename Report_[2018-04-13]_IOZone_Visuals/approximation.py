import os
import fmodpy
import numpy as np

# This directory
CWD = os.path.dirname(os.path.abspath(__file__))
CLASS_NAME = lambda obj: (repr(obj)[1:-2].split(".")[-1])
IS_NONE = lambda v: type(v) == type(None)
class BadInputShape(Exception): pass
class NoValuesProvided(Exception): pass
class NoPointsProvided(Exception): pass

# Generic class definition for creating an algorithm that can perform
# regression in the multi-variate setting. In these classes 'x' refers
# to a potentially multi-dimensional set of euclydian coordinates with
# associated single-dimensional response values referred to as 'y'.
# 
# The "fit" function is given both 'x' and associated 'y' and creates
# any necessary (or available) model for predicting 'y' at new 'x'.
# 
# The "predict" function predicts 'y' at potentially novel 'x' values.
class Approximator:
    def __init__(self, *args, **kwargs):
        # Store the extra provided keyword arguments (internal attributes)
        for name in kwargs:
            setattr(self, name, kwargs[name])

    def __str__(self):
        return CLASS_NAME(self)

    def add_point(self, x:np.ndarray, y:np.ndarray):
        pass

    def predict(self, x:np.ndarray, **kwargs):
        pass
    
    # Wrapper for 'predict' that returns a single value for a single
    # prediction, or an array of values for an array of predictions
    def __call__(self, x:np.ndarray, *args, **kwargs):
        single_response = len(x.shape) == 1
        if single_response:
            x = np.array([x])
        if len(x.shape) != 2:
            raise(BadInputShape(f"Input shape '{x.shape}' is not allowed. Only 1 or 2 dimensional inputs."))
        response = np.asarray(self.predict(x, *args, **kwargs), dtype=float)
        # Return the response values
        return response[0] if single_response else response


# Class for using the voronoi mesh to make approximations. This relies
# on some OpenMP parallelized fortran source to identify the voronoi
# mesh coefficients (more quickly than can be done in python).
class VoronoiMesh(Approximator):
    points = None
    shift = None
    scale = None
    values = None
    products = None

    # Get fortran function for calculating mesh
    def __init__(self):
        # Store the function for evaluating the mesh into self
        from voronoi import eval_mesh

    # Add a single point to the mesh
    def add_point(self, point, value):
        point = point.copy()
        # Normalize the point like others (if they have been normalized)
        if not IS_NONE(self.shift):
            point -= self.shift
        if not IS_NONE(self.scale):
            point /= self.scale
        # Update stored points
        if IS_NONE(self.points):
            self.points = np.array([point])
        else:
            self.points = np.vstack((self.points, point))
        # Update the stored values
        if IS_NONE(self.values):
            self.values = np.array(value).reshape((1,len(value)))
        else:
            self.values = np.vstack((self.values,[value]))
        # Update the fit stored in memory
        self.fit()

    # Given points, pre-calculate the inner products of all pairs.
    def fit(self, points=None, values=None, normalize=True):
        # Store values if they were provided
        if not IS_NONE(values):
            self.values = values.copy()
        # Store points if they were provided
        if not IS_NONE(points):
            self.points = points.copy()
        # Normalize points to be in unit hypercube if desired
        if normalize:
            self.shift = np.min(self.points,axis=0)
            self.points -= self.shift
            self.scale = np.max(self.points,axis=0)
            # Make sure the scale does not contain 0's
            self.scale = np.where(self.scale != 0, self.scale, 1)
            self.points /= self.scale
        # Store pairwise dot products in fortran form
        self.products = np.array(np.matmul(self.points, self.points.T), order="F")

    # Function that returns the convexified support of the
    # twice-expanded voronoi cells of fit points
    def support(self, points):
        if IS_NONE(self.points):
            raise(NoPointsProvided("Cannot compute support without points."))
        self.products = np.asarray(self.products, order="F")
        # Calculate the dot products of the points with themselves
        products = np.array(np.matmul(points, self.points.T), order="F")
        # Compute all voronoi cell values
        vvals = np.array(np.zeros((len(points), len(self.points))), order="F")
        return self.eval_mesh(self.products, products, vvals)

    # Given points to predict at, make an approximation
    def predict(self, points):
        if IS_NONE(self.values): 
            raise(NoValuesProvided("Values must be provided to make predictions."))
        # Normalize the points like others (if they have been normalized)
        if not IS_NONE(self.shift):
            points -= self.shift
        if not IS_NONE(self.scale):
            points /= self.scale
        # Calculate the support
        vvals = self.support(points)
        # Use a convex combinations of values associated with points
        # to make predictions at all provided points.
        return np.matmul(vvals, self.values)


# Wrapper class for using the Delaunay fortran code
class Delaunay(Approximator):
    points = None
    values = None
    errors = None
    # Settings for restricting the size of the interpolant
    max_points = float('inf')
    points_mass = None
    acceptable_error = 0

    # Initialize this approximator
    def __init__(self, **kwargs):
        from VTdelaunay import delaunayp
        # Store the extra provided keyword arguments (internal attributes)
        for name in kwargs:
            setattr(self, name, kwargs[name])

    # Function for adding a single point to the approximator at a time
    def add_point(self, point, value):
        # Update stored points
        if IS_NONE(self.points):
            self.points = np.array([point])
            self.values = np.array([value])
            self.points_mass = [1]
            self.errors = []
            return
        # Make a guess at the new point before adding it.
        guess = self(point)
        # Calculate the span of values
        value_shift = np.min(self.values, axis=0)
        value_scale = np.max(self.values, axis=0) - value_shift
        value_scale = np.where(value_scale != 0, value_scale, 1)
        # Calculate the normalized values
        guess_normalized = (guess - value_shift) / value_scale
        value_normalized = (value - value_shift) / value_scale
        # Calculate error (max-norm) of the normalied values
        max_norm_rel_error = max(abs(guess_normalized - value_normalized))
        # Check to see if the error is acceptable
        if (max_norm_rel_error > self.acceptable_error):
            if (len(self.points) < self.max_points):
                # Add the new known value to memory if there's room)
                self.points = np.vstack((self.points, point))
                self.values = np.vstack((self.values,[value]))
                self.points_mass += [1]
            else:
                # Otherwise, shift the nearest point towards this one
                closest_pt_idx = np.argmin((self.points - point)**2)
                mass = self.points_mass[closest_pt_idx] + 1
                # Update the point and the value associated with it
                self.points[closest_pt_idx] += (
                    point - self.points[closest_pt_idx]) / mass
                self.values[closest_pt_idx] += (
                    value - self.values[closest_pt_idx]) / mass
                # Increment the mass of this point in memory
                self.points_mass[closest_pt_idx] = mass
                # TODO:  Only shift values of vertices instead of
                #        shifting the location of the vertices.
        # For tracking purposes, record the error
        self.errors.append(max_norm_rel_error)

    # Function for immediately passing in all points and values
    def fit(self, points, values):
        self.points = None
        self.values = None
        for (pt,val) in zip(points,values):
            self.add_point(pt, val)

    # Given approximation points, return the source points and weights
    # for making a prediction at those locations
    def points_and_weights(self, approx ,interp=None, normalize=True, 
                           error_warn=True):
        # Get the interpolation points from the history
        if IS_NONE(interp): interp = self.points.copy()
        # Make sure we don't modify the approximation points
        approx = approx.copy()
        # Normalize points to be in unit hypercube if desired
        if normalize:
            shift = np.min(interp,axis=0)
            interp -= shift
            scale = np.max(interp, axis=0)
            # Make sure the scale does not contain 0's
            scale = np.where(scale != 0, scale, 1)
            interp /= scale
            # Scale the approximation points symmetrically
            approx -= shift
            approx /= scale
        # Convert interpolation and approximation points into fortran form
        interp = np.asarray(interp.T, order="F")
        approx = np.asarray(approx.T, order="F")
        # Get the predictions from VTdelaunay
        p_in = np.asarray(approx, dtype=np.float64, order="F")
        work_in = np.ones(shape=(max(5*p_in.shape[0],
                                     p_in.shape[0]*p_in.shape[1]),), 
                          dtype=np.float64, order="F")
        simp_out = np.ones(shape=(p_in.shape[0]+1, p_in.shape[1]), 
                           dtype=np.int32, order="F")
        weights_out = np.ones(shape=(p_in.shape[0]+1, p_in.shape[1]), 
                              dtype=np.float64, order="F")
        error_out = np.ones(shape=(p_in.shape[1],), 
                            dtype=np.int32, order="F")
        self.delaunayp(interp.shape[0], interp, p_in, work_in,
                       simp_out, weights_out, error_out, 
                       extrap_opt=1000)
        # Retrieve the errors from delaunay execution
        error_out = np.where(error_out != 1, error_out, 0)
        error_out = np.where(error_out != 30, error_out, 0)
        # Reset the errors to simplex of 1s (to be 0) and weights of 0s.
        bad_indices = (error_out > 1)
        num_bad = sum(bad_indices)
        if num_bad > 0:
            # Randomly assign simplices to the bad guesses
            simp_out[:,bad_indices] = np.random.randint(
                interp.shape[1], size=(simp_out.shape[0], num_bad))
            # Randomly assign convex weights to the random vertices
            weights = np.random.random(size=(num_bad, weights_out.shape[0]))
            weights /= np.sum(weights, axis=1)
            weights_out[:,bad_indices] = weights.T
        # Print error messages if desired
        if error_warn and (sum(error_out) != 0):
            unique_errors = sorted(np.unique(error_out))
            print("\n [Delaunay errors:",end="")
            for e in unique_errors:
                if (e in {0,1}): continue
                print(f" {e:3d} at {sum(error_out==e)} points", end=";")
            print("] ")
        # Adjust the output simplices and weights to be expected shape
        points  = simp_out.T - 1
        weights = weights_out.T
        return points, weights

    # Return just the points and the weights associated with them for
    # creating the correct interpolation
    def predict(self, approx, interp=None, normalize=True, error_warn=False):
        if IS_NONE(self.points):
            raise(NoPointsProvided("Points must be provided to make predictions."))
        if IS_NONE(self.values):
            raise(NoValuesProvided("Values must be provided to make predictions."))
        if (len(self.points) > len(self.points[0])):
            points, weights = self.points_and_weights(approx, interp, normalize, error_warn)
            predictions = np.array([sum(v*c for (v,c) in zip(self.values[simp],coefs))
                                    for (simp,coefs) in zip(points,weights)])
        else:
            # Use the Voronoi Mesh to evaluate under-determined structures
            model = VoronoiMesh()
            model.fit(self.points, self.values)
            predictions = model(approx)
        # Return the appropriate shaped pair of points and weights
        return predictions

# ==================================================================
#      Smoothing an approximator using independent random folds     
# ==================================================================

class FoldedApproximator(Approximator):
    k = None
    approximator = Delaunay

    # Fit multiple voronoi meshes
    def fit(self, points, values=None):
        # Default to using a k proportional to the square root num points
        k = self.k
        if IS_NONE(k): k = int(len(points) / len(points)**.5)
        # Shuffle the values and points randomly
        indices = np.arange(len(points))
        np.random.shuffle(indices)
        points = points[indices]
        if not IS_NONE(values): 
            values = values[indices]
        # Construct the individual models over subsets of points
        self.models = []
        size = int(0.5 + len(points) / k)
        for i in range(k):
            point_subset = points[i*size : (i+1)*size]
            value_subset = None
            if not IS_NONE(values):
                value_subset = values[i*size : (i+1)*size]
            self.models.append(self.approximator())
            self.models[-1].fit(point_subset, value_subset)
            
    # Predict the points by averaging the guesses over the models
    def predict(self, points):
        outputs = []
        for m in self.models:
            outputs.append( m(points) )
        return sum(outputs) / len(outputs)


# ==============================================================
#      Converting a standard approximator into a classifier     
# ==============================================================

class Classifier(Approximator):
    # Given points, pre-calculate the inner products of all pairs.
    def fit(self, points, classes=None, model=FoldedApproximator, **fit_kwargs):
        # Convert the class values into a d-dimensional vector space
        values = None
        if not IS_NONE(classes):
            self.classes = np.unique(classes)
            # Encode the d classes into a d-dimensional vector
            values = np.zeros((len(classes), len(self.classes)))
            for i,v in enumerate(self.classes):
                values[:,i] = np.where(classes == v, 1, 0)
            # Store the relative class frequency
            self.class_frequency = np.sum(values, axis=0)
            self.class_frequency /= sum(self.class_frequency)
        # Fit the approximator to the values
        self.model = model()
        self.model.fit(points, values, **fit_kwargs)

    # Use the most-likely class as the output
    def predict(self, points, frequency_compensation=False):
        # Use a convex combinations of values associated with points
        # to make predictions at all provided points.
        weighted_values = self.model.predict(points) 
        if frequency_compensation:
            weighted_values /= self.class_frequency
        class_outputs = [self.classes[i] for i in 
                         np.argmax(weighted_values, axis=1)]
        return np.array(class_outputs)


# Testing the approximators
if __name__ == "__main__":
    # Function definition
    mult = 5
    fun = lambda x: np.cos(x[0]*mult) + np.sin(x[1]*mult)
    low, upp = 0, 2
    dim = 2
    N = 500
    noise_magnitude = 1
    np.random.seed(0)

    # Settings for testing
    random = True
    add_noise = True
    classification = False
    # surf = Delaunay()
    surf = VoronoiMesh()
    # surf = FoldedApproximator(k=2, approximator=Delaunay)
    # surf = Classifier()

    # Generate the x-points 
    if random:
        x = np.random.random(size=(N,dim)) * (upp-low) + low
    else:
        N = int(round(N ** (1/dim)))
        x = np.array([r.flatten() for r in np.meshgrid(np.linspace(low,upp,N), np.linspace(low,upp,N))]).T

    # Calculate the function response at all points
    y = np.array([fun(v) for v in x])
    # Add noise to the function evaluation
    if add_noise:
        y += np.random.random(size=y.shape) * noise_magnitude
    # Convert to a classification problem if necessary
    if classification:
        small = min(y)
        large = max(y)
        y = np.where(y < 1/3*(large-small)+small, 1, 0) + \
            (np.where((1/3*(large-small)+small <= y), 2,0)*np.where(y < 2/3*(large-small)+small,1,0)) + \
            np.where(2/3*(large-small)+small <= y, 3, 0)

    print("Fitting the approximator...")
    surf.fit(x,y)
    print("  done")
    from util.plotly import Plot
    print("Plotting the approximation...")
    # Plot the training points and the surface
    p = Plot()
    p.add("Training Points", *x.T, y)
    p.add_func(str(surf), surf, *([(low,upp)]*dim), plot_points=2000,
               use_gradient=True)
    print("Creating HTML...")
    p.plot(z_range=[-3.1,3.1], file_name="test_plot.html")
    print("  done...")

