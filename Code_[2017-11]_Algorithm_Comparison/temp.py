from util.system import load
from util.data import Data


name = "del_error_61.pkl"
# name = "delaunay_failure_{unique_errors[0]}.pkl"
params = load(name)

pts_shape_0, pts_shape_1_, pts_in, p_in_shape_1, p_in, simp_out, weights_out, error_out, eps, pmode, chunksize = params


print("Train points shape:", pts_in.shape)
print("Test points shape: ", p_in.shape)

from util.algorithms import DelaunayPNC10 as Del

for i in range(10000):
    print(i,end="\r")
    model = Del(chunksize=5)
    model.fit(pts_in.T)
    model(p_in.T)


# train = Data()
# for p in pts_in.T:
#     train.append(p)

# test = Data()
# for p in p_in.T:
#     test.append(p)

# train.save("err_61_at_97-98_train.csv")
# test.save("err_61_at_97-98_test.csv")
 
