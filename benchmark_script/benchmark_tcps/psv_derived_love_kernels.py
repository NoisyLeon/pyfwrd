import sys
sys.path.append('/home/leon/code/pysurf')

"""
Conclusion:

Flat Earth:
1. dh = 1., 2., 4.; derived kernels are correct.
    But the larger the dh, the larger the difference between the directly-computed and derived kernels!

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
# tcps2 = tcps.tcps_solver(m)
# tcps2.init_default(nl=100, dh=2.)
# tcps2.verbose=1
# tcps2.solve_PSV()
# 
# tcps3 = tcps.tcps_solver(m)
# tcps3.init_default(nl=50., dh=4.)
# tcps3.verbose=1
# tcps3.solve_PSV()

# 
# 
i=5
# # plt.plot(tcps1.dArr.cumsum(), np.abs(tcps1.dcdah[i, :]-tcps1.dcdah1[i, :])/tcps1.dcdah[i, :], 'ro-', ms=10)
# # plt.plot(tcps1.dArr.cumsum(), np.abs(tcps1.dcdah[i, :]-tcps1.dcdah2[i, :])/tcps1.dcdah[i, :], 'bo-', ms=10)
#
# for i in xrange(10):
#     plt.plot(tcps1.dArr.cumsum(), tcps1.dcdA[i, :], 'o-', ms=5)
#     plt.plot(tcps1.dArr.cumsum(), tcps1.dcdA1[i, :], '^-', ms=5)
# plt.plot(tcps1.dArr.cumsum(), tcps1.dcdA[i, :], 'o-', ms=5)
# plt.plot(tcps1.dArr.cumsum(), tcps1.dcdA1[i, :], '^-', ms=5)

plt.plot(tcps1.dArr.cumsum(), tcps1.dcdC[i, :], 'o-', ms=5)
plt.plot(tcps1.dArr.cumsum(), tcps1.dcdC1[i, :], '^-', ms=5)
# plt.plot(tcps2.dArr.cumsum(), tcps2.dcdbv[i, :], 'o-', ms=10)
# plt.plot(tcps2.dArr.cumsum(), tcps2.dcdbv1[i, :], '^-', ms=10)
# plt.plot(tcps3.dArr.cumsum(), tcps3.dcdbv[i, :], 'o-', ms=10)
# plt.plot(tcps3.dArr.cumsum(), tcps3.dcdbv1[i, :], '^-', ms=10)
plt.show()

