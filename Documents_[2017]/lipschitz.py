import numpy as np

# Class for creating a surface with minimal lipschitz error
class Lipschitz:
    def __init__(self):
        self.constant = 0
        self.nodes = []
        self.values = []

    # Add a set of points to this class for storage
    def add_points(self, points):
        for pt in points:
            # Update the stored lipschitz constant
            for (n,v) in zip(self.nodes,self.values):
                dist = np.sqrt(np.sum((n-pt[:-1])**2))
                self.constant = max(self.constant, abs(pt[-1]-v)/dist)
            # Add the point to the set of nodes
            self.nodes.append(pt[:-1])
            self.values.append(pt[-1])

    # Return the value with equal plus and minus lipschitz error
    def __call__(self, x):
        # Initialize holder for plus and minus error
        low_upp = [-float('inf'),float('inf')]
        for (n,v) in zip(self.nodes, self.values):
            dist = np.sqrt(np.sum( (n - x)**2 ))
            low_upp[0] = max(low_upp[0], v - dist*self.constant)
            low_upp[1] = min(low_upp[1], v + dist*self.constant)
        # Remove all infinities from surface
        while float('inf') in low_upp:
            low_upp[low_upp.index(float('inf'))] = 0
        # Return the average of the lower and upper values
        return low_upp

