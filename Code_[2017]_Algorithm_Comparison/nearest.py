import random
import numpy as np
from scipy.optimize import minimize

SMALL_NUM = .00001
BIG_NUM = 1000000

L2 = lambda xs,x: np.sqrt(np.sum((xs-x)**2,axis=1))
L2SQ = lambda xs,x: np.sum((xs-x)**2, axis=1)

LINEAR = lambda r: max(0,(1+SMALL_NUM) - r)
QUADRATIC = lambda r: max(0,(1+SMALL_NUM) - r)**2

class Linn:
    # Initialize the Linn regressor.
    def __init__(self, radius=LINEAR, distance=L2):
        self.radius = radius
        self.distance = distance

    # Fit a set of points
    def fit(self, points, values):
        self.points = points.copy()
        self.values = values.copy()

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

    # Give the weights associated with a particular input point
    def weights(self, xs, nearest_ind=None):
        if type(nearest_ind) == type(None):
            nearest_ind = random.randint(0,len(self.points)-1)
        predictions = []
        for x in xs:
            nearest = self.points[nearest_ind]
            dist_to_point = self.distance([nearest], x)[0]
            # Get the unit-length vector from (nearest -> x)
            vector = (x - nearest)
            vector /= np.sum(vector**2)**.5

            #      Identify boundary     
            # ===========================
            mult = 1
            shift = 1
            while (SMALL_NUM < shift) and (mult < BIG_NUM):
                neighbor_ind = np.argmin(self.distance(self.points, nearest+mult*vector))
                if (neighbor_ind == nearest_ind):
                    shift *= 2
                else:
                    mult -= shift
                    shift /= 2
                mult += shift
            # NOTE: Intentionally added an extra 'shift' after
            #       stopping criteria to ensure we are past boundary.

            # Calculate the distance from the nearest to the boundary
            dist_to_boundary = self.distance([nearest], nearest+mult*vector)[0]

            # Get the maximum distance from center with non-zero weight
            max_dist = 2*dist_to_boundary # self.distance([nearest], nearest+mult*vector)[0]
            # Calculate the weight and record it
            if mult >= BIG_NUM:
                weight = 1.0
            else:
                weight = self.radius(dist_to_point / max_dist)
            predictions.append(weight)

        return predictions

        

    # Generate a prediction for a new point
    def predict(self, xs):
        predictions = []
        for x in xs:
            # nearest_ind = 0 # np.argmin(self.distance(self.points, x))
            weights = []
            for nearest_ind in range(len(self.points)):
                nearest = self.points[nearest_ind]
                dist_to_point = self.distance([nearest], x)[0]
                # Get the unit-length vector from (nearest -> x)
                vector = (x - nearest)
                vector /= np.sum(vector**2)**.5

                #      Identify boundary     
                # ===========================
                mult = 1
                shift = 1
                while (SMALL_NUM < shift) and (mult < BIG_NUM):
                    neighbor_ind = np.argmin(self.distance(self.points, nearest+mult*vector))
                    if (neighbor_ind == nearest_ind):
                        shift *= 2
                    else:
                        mult -= shift
                        shift /= 2
                    mult += shift
                # NOTE: Intentionally added an extra 'shift' after
                #       stopping criteria to ensure we are past boundary.

                # Calculate the distance from the nearest to the boundary
                dist_to_boundary = self.distance([nearest], nearest+mult*vector)[0]

                # Get the maximum distance from center with non-zero weight
                max_dist = 2*dist_to_boundary # self.distance([nearest], nearest+mult*vector)[0]
                # Calculate the weight and record it
                if mult >= BIG_NUM:
                    weights.append( 1.0 )
                else:
                    weights.append( self.radius(dist_to_point / max_dist) )

            predictions.append( sum(np.array(weights)*self.values)/sum(weights) )
            # predictions.append(weights[0])

        return predictions

