import sys
sys.path.append('/home/leon/code/pysurf')

"""
Conclusion:


1. Partial derivatives, dut/dz have been benchmarked!
2. eigenfunctions, as expected, do not change with layer thickness !
3. group velocity computed using (k*I1)/omege/I0 is NOT very accurate for long periods

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
tcps1.verbose=1
tcps1.solve_SH()

# 
tcps2 = tcps.tcps_solver(m)
tcps2.init_default(nl=100., dh=2.)
tcps2.verbose=1
tcps2.solve_SH()
# 
# 
# i=5
# plt.plot(tcps1.dArr.cumsum(), tcps1.dutdz[i, :], 'ro-', ms=10)
# plt.plot(tcps1.dArr[1:].cumsum(), tcps1.dutdz1[i, :], 'bo-', ms=10)
# #     
# for i in xrange(10):
#     # plt.figure()
#     plt.plot(tcps1.dArr.cumsum(), tcps1.ut[i, :], 'ro-', ms=10)
#     plt.plot(tcps2.dArr.cumsum(), tcps2.ut[i, :], 'bo-', ms=10)
    # plt.plot(tcps1.dArr.cumsum(), tcps1.dutdz[i, :], 'ro-', ms=10)
    # plt.plot(tcps1.dArr[1:].cumsum(), tcps1.dutdz1[i, :], 'bo-', ms=10)
# # # # 
# # # # plt.figure()
plt.plot(tcps1.T, tcps1.Vgr, 'ro', ms=10)
plt.plot(tcps1.T, tcps1.Vgr1, 'b^', ms=10)
plt.show()

