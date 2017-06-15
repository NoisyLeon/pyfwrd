import eigen
import vmodel
import numpy as np

m=vmodel.model1d()
m.get_radius(400., 1.)
m.model_prem()
eig = eigen.eigen_solver(m)
eig.initialization()
