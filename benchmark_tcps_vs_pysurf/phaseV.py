import sys
sys.path.append('/home/leon/code/pysurf')

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

tcps2 = tcps.tcps_solver(m)
tcps2.init_default_2()
tcps2.verbose=1
tcps2.solve_PSV()


# plt.plot((eig1.T), (eig1.Vph[0,:]/1000.), 'ro', ms=10)
# plt.plot(tcps1.T, tcps1.Vph, 'bx', ms=10)
# plt.show()
