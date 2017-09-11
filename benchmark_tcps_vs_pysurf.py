import eigen, tcps
import vmodel
import numpy as np
import matplotlib.pyplot as plt

m=vmodel.model1d()
m=vmodel.model1d()
m=vmodel.read_model(m, 'ak135.txt')
# m.earth_flattening()
eig1 = eigen.eigen_solver(m)
eig1.init_default()
eig1.solve_PSV()


tcps1 = tcps.tcps_solver(m)
tcps1.init_default()
tcps1.solve_PSV()


plt.plot((eig1.T), (eig1.Vph[0,:]/1000.), 'o')
plt.plot(tcps1.T, tcps1.Vph, 'x')
plt.show()
