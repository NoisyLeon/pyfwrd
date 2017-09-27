import sys
sys.path.append('/home/leon/code/pysurf')

"""
Conclusion:

Flat Earth:
1. T = 5 - 100. sec, dT = 5 sec
    Results are perfectly consistent, except for the fact T = 5 sec, aniprop yield results at T = 10 sec.

"""
import eigen, tcps, aniproppy
import vmodel
import numpy as np
import matplotlib.pyplot as plt

m=vmodel.model1d()
m.model_ak135_cps()
m.flat=1

tcpsR = tcps.tcps_solver(m)
tcpsR.init_default()
tcpsR.solve_PSV()


tcpsL = tcps.tcps_solver(m)
tcpsL.init_default()
tcpsL.solve_SH()
# 
ani  = aniproppy.aniprop_solver(m)
ani.init_default_2()
ani.solve()
# # 
# tcps3 = tcps.tcps_solver(m)
# tcps3.init_default(nl=50., dh=4.)
# tcps3.verbose=1
# tcps3.solve_PSV()
# 
# # 
# # 
# i=5
# # # plt.plot(tcps1.dArr.cumsum(), np.abs(tcps1.dcdah[i, :]-tcps1.dcdah1[i, :])/tcps1.dcdah[i, :], 'ro-', ms=10)
# # # plt.plot(tcps1.dArr.cumsum(), np.abs(tcps1.dcdah[i, :]-tcps1.dcdah2[i, :])/tcps1.dcdah[i, :], 'bo-', ms=10)
# # 
# plt.plot(tcps1.dArr.cumsum(), tcps1.dcdbv[i, :], 'o-', ms=10)
# plt.plot(tcps1.dArr.cumsum(), tcps1.dcdbv1[i, :], '^-', ms=10)
# plt.plot(tcps2.dArr.cumsum(), tcps2.dcdbv[i, :], 'o-', ms=10)
# plt.plot(tcps2.dArr.cumsum(), tcps2.dcdbv1[i, :], '^-', ms=10)
# plt.plot(tcps3.dArr.cumsum(), tcps3.dcdbv[i, :], 'o-', ms=10)
# plt.plot(tcps3.dArr.cumsum(), tcps3.dcdbv1[i, :], '^-', ms=10)
# plt.show()