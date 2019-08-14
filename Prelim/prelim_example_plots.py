from util.approximate.testing import test_plot

from util.approximate import MARS, SVR, NeuralNetwork, \
    Delaunay, ShepMod, LSHEP


# model = MARS
# m = model(max_bases=5, max_interaction=1)
# p, x, y = test_plot(m, random=True, N=20)
# p.show(file_name="demo_MARS.html")

# model = SVR
# m = model()
# p, x, y = test_plot(m, random=True, N=200)
# p.add("20 train", *(x[:20].T), y[:20], color=p.color(0))
# p.show(file_name="demo_SVR.html")

# model = NeuralNetwork
# m = model(layers=(5,3))
# p, x, y = test_plot(m, random=True, N=20)
# p.show(file_name="demo_MLP.html")

# model = Delaunay
# m = model()
# p, x, y = test_plot(m, random=True, N=20)
# p.show(file_name="demo_Delaunay.html")

model = ShepMod
m = model()
p, x, y = test_plot(m, random=True, N=20)
p.show(file_name="demo_ShepMod.html")

# model = LSHEP
# m = model(radius=4)
# p, x, y = test_plot(m, random=True, N=20)
# p.show(file_name="demo_LSHEP.html")
