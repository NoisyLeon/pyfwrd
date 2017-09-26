"""
The kernels are slightly different!
"""

import tcps
import vmodel
import numpy as np
import matplotlib.pyplot as plt

m=vmodel.model1d()
m.model_ak135_cps()
m.flat=1

tcps1 = tcps.tcps_solver(m)
tcps1.init_default()
tcps1.solve_PSV()

m=vmodel.model1d()
m.model_ak135_cps()
m.earth_flattening()
tcps2 = tcps.tcps_solver(m)
tcps2.init_default()
tcps2.solve_PSV()


i=5
plt.plot(tcps1.dArr.cumsum(), tcps1.dcdah[i, :], 'ro-', ms=10)
plt.plot(tcps2.dsph.cumsum(), tcps2.dcdah[i, :], 'bo-', ms=10)
plt.show()
