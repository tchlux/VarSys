# Author and Email: 
#   Thomas Lux (tchlux@vt.edu)
# 
# Modifications:
#   2018 April -(TL)- Created 'regular_simplex' function.
# 
# Given "d" categories that need to be converted into real space,
# generate a regular simplex in (d-1)-dimensional space. (all points
# are equally spaced from each other and the origin) This process
# guarantees that all points are not placed on a sub-dimensional
# manifold (as opposed to "one-hot" encoding).
def regular_simplex(num_categories):
    import numpy as np
    class InvalidNumberOfCategories(Exception): pass
    # Special cases for one and two categories
    if num_categories <= 1:
        raise(InvalidNumberOfCategories(
            "Number of categories must be an integer greater than 1."))
    elif num_categories == 2:
        return np.array([[0.],[1.]])
    # Standard case for >2 categories
    d = num_categories
    # Initialize all points to be zeros.
    points = np.zeros((d,d-1))
    # Set the initial first point as 1
    points[0,0] = 1
    # Calculate all the intermediate points
    for i in range(1,d-1):
        # Set all points to be flipped from previous calculation while
        # maintaining the angle "arcos(-1/d)" between vectors.
        points[i:,i-1] = -1/(d-i) * points[i-1,i-1]
        # Compute the new coordinate using pythagorean theorem
        points[i,i] = (1 - sum(points[i,:i]**2))**(1/2)
    # Set the last coordinate of the last point as the negation of the previous
    points[i+1,i] = -points[i,i]
    # Return the regular simplex
    return points
