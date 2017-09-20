import sys
sys.path.append('/home/leon/code/pysurf')

"""
Conclusion:

Flat Earth:
1. 
    

Spherical Earth:
1. dh = 1., 2., 4.; derived kernels are correct
    But the larger the dh, the larger the difference between the directly-computed and derived kernels!

"""
import eigen, tcps
import vmodel
import numpy as np
import matplotlib.pyplot as plt

m=vmodel.model1d()
m.model_ak135_cps()
m.flat=1

tcps1 = tcps.tcps_solver(m)
tcps1.init_default()
tcps1.solve_PSV()

# 
tcps2 = tcps.tcps_solver(m)
m.add_perturb_layer(0., 10., 4, -0.6,True)
tcps2.init_default()
# tcps2.verbose=1
tcps2.solve_PSV()

tcps1.psv_perturb_disp_vel(tcps2)
# tcps3 = tcps.tcps_solver(m)
# tcps3.init_default(nl=50., dh=4.)
# tcps3.verbose=1
# tcps3.solve_PSV()



plt.plot(tcps1.T, tcps1.Vph, 'o', ms=10)
plt.plot(tcps1.T, tcps1.Vph_pre, '^', ms=10)
plt.plot(tcps2.T, tcps2.Vph, 'x', ms=10)

plt.show()

