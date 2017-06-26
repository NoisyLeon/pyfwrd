import eigen
import vmodel
import numpy as np
import matplotlib.pyplot as plt

m=vmodel.model1d()
# m.get_radius(4000., 1.)
# m.model_prem()
m=vmodel.model1d()
m=vmodel.read_model(m, 'ak135.txt')
eig1 = eigen.eigen_solver(m)
eig1.initialization()
# eig1.solve_PSV()

# eig2 = eigen.eigen_solver(m)
# eig2.init2()
# eig2.solve_PSV()
# 
# plt.plot(eig1.r, eig1.r2data[0, 8, :], 'o', ms=5)
# plt.plot(eig2.r, eig2.r2data[0, 8, :], '^', ms=5)
# plt.show()
