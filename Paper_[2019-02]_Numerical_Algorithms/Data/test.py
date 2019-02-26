from util.math import SMALL
from util.approximate import Delaunay
import numpy as np

print()
print("Loading hard problem..")

a = np.loadtxt("hard_problem.csv", delimiter=",", skiprows=1)
train, test = a[:-1,:], a[-1,:]

train = np.asfortranarray(train.T)
test = np.asfortranarray(test[:,None])

pts_in = train
p_in = test
simp_out = np.ones(shape=(p_in.shape[0]+1, p_in.shape[1]), 
                   dtype=np.int32, order="F")
weights_out = np.ones(shape=(p_in.shape[0]+1, p_in.shape[1]), 
                      dtype=np.float64, order="F")
error_out = np.ones(shape=(p_in.shape[1],), 
                    dtype=np.int32, order="F")

print("Calling delaunay serial subroutine manually..")

m = Delaunay()
m.delaunays(train.shape[0], train.shape[1],
            pts_in, p_in.shape[1], p_in, simp_out,
            weights_out, error_out, extrap=100.0, ibudget=1000,
            eps=2*SMALL)

print("done.")
