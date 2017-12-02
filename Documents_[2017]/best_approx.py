import numpy as np
from plotly_interface import *
from scipy import integrate

# Fint the hyper-plane that splits the data such that the CDF's of the
# response on each side of the hyper-plane have the largest area
# between the curves. That means that we are maximizing:
#    integral( (cdf_{left} - cdf_{right})^2 )



class Approximator:
    # Assumes that rows of the data are points to interpolated, and
    # that the last column is the value that should be interpolated
    def __init__(self, data):
        self.data = data
        self.weights = np.ones(data.shape[0])
        self.model = [np.ones(data.shape[1]-1)]
        self.min_val = min(self.data[:,-1])
        self.max_val = min(self.data[:,-1])
        
    # Return the negative 
    def evaluator(self, plane):
        side = np.where(np.dot(self.data[:,:-1], plane) > 0, 1, 0)    
        cdf1 = cdf_fit_func(self.data[np.where(side == 1)])
        cdf2 = cdf_fit_func(self.data[np.where(side == 0)])    
        return -quad(lambda x: abs(cdf1(x) - cdf2(x)), self.min_val, self.max_val)[0]


    def improve(self):
        minimize(

    def __eval__(self, x):
        prod = 1.0

        

if __name__ == __main__:

    # Testing code
