import numpy as np
from plotly_interface import Plot

data = np.array([
    [1.2, 1.2],
    [3.6, 3.1],
    [4.3,-1.3],
    [6.1, 2.7],
    [7.8, 1.4]
])

def u1(x):
    return abs(x)

def u2(x):
    return np.exp(-x**2)

def F(x,u=u1,d=data):
    return sum(cj*u(x - xj) for (xj,cj) in d)

p = Plot()

p.add("Data", data[:,0], data[:,1])
p.add_func("Abs", lambda x: F(x,u1), [min(data[:,0])-1, max(data[:,0])+1])
p.add_func("e", lambda x: F(x,u2), [min(data[:,0])-1, max(data[:,0])+1])
p.plot(fixed=False)
